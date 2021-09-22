
#Note:
#This function (plot_PFS_Microbiome_curve) has to be used with the JAMS package (https://github.com/johnmcculloch/JAMS_BW). It has been used for plotting a curve showing the separations of P and NP groups based on Progression Free Survival throughout time.

plot_PFS_Microbiome_curve <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, algorithm = "PCA", distmethod = "bray", colourby = NULL, shapeby = NULL, sizeby = NULL, pairby = NULL, dotsize = 2, dotborder = NULL, log2tran = TRUE, transp = TRUE, perplx = NULL, max_neighbors = 15, ellipse = FALSE, plotcentroids = FALSE, show_centroid_distances = FALSE, addtit = NULL, cdict = NULL, grid = TRUE, forceaspectratio = NULL, threads = 1, class_to_ignore = "N_A", Slice_thickness_days = 30, convert_days_to_months = FALSE, crop_to_time = NULL, return_Ordinations = TRUE, return_metadata_slices = FALSE, return_Slice_df = FALSE, stratifybytreatment = TRUE, loess_span = 0.25, permanova_permutations = 10000, add_OR_asterisk = TRUE, NP_FollowUp_column = "PFS_days"){


    set.seed(4140)

    subsetby = NULL

    #Remove samples bearing categories within class_to_ignore
    #valid_vars <- c(colourby, shapeby, sizeby, "FollowUp_Time_days")[which(!is.na(c(colourby, shapeby, sizeby, "FollowUp_Time_days")))]
    valid_vars <- c(colourby, shapeby, sizeby)[which(!is.na(c(colourby, shapeby, sizeby)))]

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = valid_vars, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

    subset_points <- "none"
    sp = 1

    #Create list vector to hold plots
    gvec <- NULL
    gvec <- list()
    pcatitbase <- paste(algorithm, "of", analysisname)

    samplesToKeep <- rownames(colData(obj))
    subsetname <- "no_sub"
    pcatit <- pcatitbase

    pcatit <- paste(c(pcatit, presetlist$filtermsg), collapse = "\n")

    currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

    currpt <- as.data.frame(colData(currobj))

    relevantcontinuousvars <- c("Days_on_treatment", "FollowUp_Time_days", "Progression_Event", "PFS_days", "Study_Clin_Response", "ORR")[c("Days_on_treatment", "FollowUp_Time_days", "Progression_Event", "PFS_days", "Study_Clin_Response", "ORR") %in% colnames(currpt)]

    currpt <- currpt[ , c("Sample", valid_vars, relevantcontinuousvars)]

    #Regenerate Slices
    numslices <- as.integer(max(as.numeric(currpt$FollowUp_Time_days[which(currpt$FollowUp_Time_days != "N_A")])) / Slice_thickness_days)

    Slices <- (1:numslices) * Slice_thickness_days

    #Add exactly 1, 2, and 3 years
    #anniversaries <- c(1:3 * 365)[c(1:3 * 365) <= max(as.numeric(currpt$FollowUp_Time_days))]
    #Slices <- sort(unique(c(anniversaries, Slices)))

    if (all(c(stratifybytreatment, ("Days_on_treatment" %in% relevantcontinuousvars)))){
        stratifybytreatment <- TRUE
        progressioncats <- c("Under_Treatment_No_Progression", "Post_Treatment_No_Progression", "Progression")
    } else {
        stratifybytreatment <- FALSE
        progressioncats <- c("No_Progression", "Progression")
    }

    flog.warn(paste("FollowUp data for Non-Progressors will be picked up from the", NP_FollowUp_column, "column"))

    for (days in Slices){
        Var_Name <- paste("PFS", days, sep = "_")
        #Start by censoring everyone and then fill each one with appropriate response at that time
        currpt[ , Var_Name] <- NA
        #Loop through each sample
        for (spos in 1:length(currpt[ , Var_Name])){
            #If progressed, see if progression had occurred at that timepoint
            if (as.numeric(currpt$Progression_Event[spos]) == 1){
                if (as.numeric(currpt$PFS_days[spos]) < days){
                    currpt[spos, Var_Name] <- "Progression"
                } else {
                    if (stratifybytreatment){
                        if (as.numeric(currpt$Days_on_treatment[spos]) >= days){
                            currpt[spos, Var_Name] <- "Under_Treatment_No_Progression"
                       } else {
                            currpt[spos, Var_Name] <- "Post_Treatment_No_Progression"
                        }
                    } else {
                        currpt[spos, Var_Name] <- "No_Progression"
                    }
                }

            } else {
                #There is no progression, so state No_Progression if followup days is longer than current time slice. If it is shorter, then leave censored.

                if (as.numeric(currpt[spos, NP_FollowUp_column]) >= days){
                    if (stratifybytreatment){
                        if (as.numeric(currpt$Days_on_treatment[spos]) >= days){
                            currpt[spos, Var_Name] <- "Under_Treatment_No_Progression"
                        } else {
                            currpt[spos, Var_Name] <- "Post_Treatment_No_Progression"
                        }
                    } else {
                        currpt[spos, Var_Name] <- "No_Progression"
                    }
                }
            }
        }
    }

    #Fix NAs to JAMS "N_A"
    for (colm in 1:ncol(currpt)){
        #Coerce to character
        currpt[ , colm] <- as.character(currpt[ , colm])
        currpt[is.na(currpt[ , colm]) , colm] <- "N_A"
    }

    #eliminate Slices with low proportions of response
    Slicenames <- paste("PFS", Slices, sep = "_")
    Slicedf <- data.frame()
    for (slc in Slicenames){
        Slicedf[slc, 1:length(c("N_A", progressioncats))] <- (as.numeric(table(currpt[ , slc])[c("N_A", progressioncats)]))
    }

    colnames(Slicedf) <- c("Censored", progressioncats)
    Slicedf$Censored[is.na(Slicedf$Censored)] <- 0

    if (stratifybytreatment){
        Slicedf$ProgressionRatio <- round((Slicedf$Progression / rowSums(Slicedf[,1:4])), 2)
    } else {
        Slicedf$ProgressionRatio <- round((Slicedf$Progression / rowSums(Slicedf[,1:3])), 2)
    }

    Slicedf$CensoringRatio <- round((Slicedf$Censored / nrow(currpt)), 2)
    Slicedf$Days <- Slices

    #Eliminate rows with NA
    Slicedf <- Slicedf[!is.na(Slicedf$ProgressionRatio), ]
    #Eliminate proportions of Progression smaller than 5%
    Slicedf <- subset(Slicedf, ProgressionRatio >= 0.05)
    Slicedf <- subset(Slicedf, ProgressionRatio <= 0.95)
    Slicedf <- subset(Slicedf, CensoringRatio <= 0.95)
    Slicedf$Months <- round((Slicedf$Days / 30), 1)

    #Add number of patients in treatment at each slice
    #if ("Days_on_treatment" %in% colnames(currpt)){
    #   Slicedf$Patients_under_treatment <- sapply(Slicedf$Days, function (x) { length(which(as.numeric(currpt$Days_on_treatment) >= x)) })
    #}

    ValidSlices <- rownames(Slicedf)
    currpt <- currpt[ , c("Sample", valid_vars, relevantcontinuousvars, ValidSlices)[c("Sample", valid_vars, relevantcontinuousvars, ValidSlices) %in% colnames(currpt)]]

    #Get counts matrix
    countmat <- as.matrix(assays(currobj)$BaseCounts)

    #Protect against rows with empty data
    rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
    countmat <- countmat[rowsToKeep, ]

    if (ignoreunclassified == TRUE){
        dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
        rowsToKeep <- which(!(rownames(countmat) %in% dunno))
        countmat <- countmat[rowsToKeep, ]
    }

    #log2 transform if applicable
    if (log2tran == TRUE){
        #Transform to log2 space
        countmat <- convert_matrix_log2(mat = countmat, transformation = "to_log2")
        pcatit <- paste(c(pcatit, "Matrix log2 transformed"), collapse = "\n")
    }

    n <- nrow(countmat)
    comp <- 1:60
    rowVars <- rowSds(countmat)
    countmat <- countmat[order(rowVars, decreasing = TRUE), ]

    countmat <- t(countmat)

    ##########
    ## Get permanova for each slice

    Slicedf$Permanova <- NA
    Slicedf$SumsOfSqs <- NA
    Slicedf$MeanSqs <- NA
    Slicedf$FModel <- NA
    Slicedf$R2 <- NA
    Slicedf$Median_dist <- NA
    Slicedf$Mean_dist <- NA

    for (slc in ValidSlices){
        set.seed(4140)
        #Get countmat with uncensored samples
        ValidSamples <- currpt$Sample[(currpt[ , slc] != "N_A")]
        slc_countmat <- countmat[ValidSamples, ]
        d <- vegdist(slc_countmat, method = distmethod, na.rm = TRUE)
        slc_pheno <- currpt[ValidSamples, ]
        if (stratifybytreatment){
            slc_pheno[which(slc_pheno[ , slc] == "Post_Treatment_No_Progression"), slc] <- "No_Progression"
            slc_pheno[which(slc_pheno[ , slc] == "Under_Treatment_No_Progression"), slc] <- "No_Progression"
        }

        Slicedf[slc, "Permanova"] <- vegan::adonis(as.formula(paste("d ~ ", slc)), data = slc_pheno, permutations = permanova_permutations)$aov.tab$`Pr(>F)`[1]
        Slicedf[slc, "SumsOfSqs"] <- vegan::adonis(as.formula(paste("d ~ ", slc)), data = slc_pheno, permutations = permanova_permutations)$aov.tab$SumsOfSqs[1]
        Slicedf[slc, "MeanSqs"] <- vegan::adonis(as.formula(paste("d ~ ", slc)), data = slc_pheno, permutations = permanova_permutations)$aov.tab$MeanSqs[1]
        Slicedf[slc, "FModel"] <- vegan::adonis(as.formula(paste("d ~ ", slc)), data = slc_pheno, permutations = permanova_permutations)$aov.tab$`F.Model`[1]
        Slicedf[slc, "R2"] <- vegan::adonis(as.formula(paste("d ~ ", slc)), data = slc_pheno, permutations = permanova_permutations)$aov.tab$R2[1]
        #ANOSIM
        Slicedf[slc, "ANOSIM"] <- anosim(slc_countmat, as.factor(slc_pheno[ , slc]), distance = distmethod, permutations = permanova_permutations)$signif


        #Get straight median dissimilarities, without reducing ordination
        curr_progressors <- rownames(slc_pheno)[which(slc_pheno[ , slc] == "Progression")]
        curr_no_progressors <- rownames(slc_pheno)[which(slc_pheno[ , slc] == "No_Progression")]

        Slicedf[slc, "Median_dist"] <- median(as.matrix(d)[curr_progressors, curr_no_progressors])
        Slicedf[slc, "Mean_dist"] <- mean(as.matrix(d)[curr_progressors, curr_no_progressors])
        #Slicedf[slc, "Min_dist"] <- min(as.matrix(d)[curr_progressors, curr_no_progressors])
    }

    Slicedf$Reciprocal_Permanova <- round((1 / Slicedf$Permanova), 3)
    Slicedf$Reciprocal_ANOSIM <- round((1 / Slicedf$ANOSIM), 3)

    Slicedf$Euclidean_Distance_P_NP_2D <- NA
    Slicedf$Euclidean_Distance_P_NP_OmniD <- NA

    #Get master ordination
    if (algorithm == "tUMAP"){
        set.seed(1234)
        n_neighbors <- min((nrow(countmat) - 1), max_neighbors)
        tumap_out <- tumap(countmat, n_components = 2, n_neighbors = n_neighbors, verbose = FALSE, n_threads = threads, init = "spca")
        dford <- as.data.frame(tumap_out)
        rownames(dford) <- rownames(currpt)
        colnames(dford) <- c("PC1", "PC2")
        xl <- "tUMAP 1"
        yl <- "tUMAP 2"

    } else {
        set.seed(1234)
        d <- vegdist(countmat, method = distmethod, na.rm = TRUE)
        pcaRes <- prcomp(d)
        ord <- pcaRes$x
        comps_have <- colnames(ord)
        vars <- pcaRes$sdev^2
        vars <- round(vars/sum(vars), 5) * 100

        xl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[1]], vars[comp[1]])
        yl <- sprintf("%s: %.2f%% variance", colnames(ord)[comp[2]], vars[comp[2]])

        dford <- as.data.frame(ord)

        #Make a data frame with how variance is explained by which components
        vardf <- data.frame(Component = colnames(ord), Variance = vars, Cumulative_variance = cumsum(vars))
        vardf$Component <- factor(vardf$Component, levels = vardf$Component)

        varplot <- ggplot(vardf, aes(x = Component)) + geom_bar(aes(y = Variance), fill = 'blue', stat = "identity") + geom_point(aes(y = Cumulative_variance), colour = "black", pch = 16, size = 1) + geom_path(aes(y = Cumulative_variance, group = 1)) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.6, size = rel(0.5))) + labs(title = "Variance explained by each PCoA component", x = 'Component', y = 'Variance')

        #print(varplot)

    }

    dford$Sample <- rownames(dford)
    dford <- left_join(dford, currpt, by = "Sample")
    dford$FollowUp_Time_days <- as.numeric(dford$FollowUp_Time_days)
    rownames(dford) <- dford$Sample

    #Make a list of plots for each ValidSlice
    SliceOrdinations <- list()
    dford_plot <- dford
    if (stratifybytreatment){
        for (colm in 1:ncol(dford_plot)){
            dford_plot[which(dford_plot[ , colm] == "Post_Treatment_No_Progression"), colm] <- "No_Progression"
            dford_plot[which(dford_plot[ , colm] == "Under_Treatment_No_Progression"), colm] <- "No_Progression"
        }
    }

    #Get distances between centroids
    for (slc in ValidSlices){
        slc_dford <- dford_plot
        slc_dford[which(slc_dford[ , slc] == "N_A"), slc] <- "Censored"
        colnames(slc_dford)[which(colnames(slc_dford) == slc)] <- "Colours"

        centroids <- aggregate(cbind(PC1, PC2) ~ Colours, slc_dford, mean)
        colnames(centroids)[c(2, 3)] <- c("meanPC1", "meanPC2")

        centroiddf <- left_join(slc_dford, centroids, by = "Colours")
        rownames(centroids) <- centroids$Colours

        #Get Euclidean distance between centroids
        centroiddist <- as.matrix(dist(centroids[,2:ncol(centroids)], method = "euclidean"))

        #Pump value into Slicedf
        Slicedf[slc , "Euclidean_Distance_P_NP_2D"] <- centroiddist["Progression", "No_Progression"]

        if (all(paste0("PC", 1:3) %in% colnames(slc_dford))){
            cat(paste("Calculating omni-dimensional aggregates in", slc, "\n"))

            centroidsOmniD <- aggregate(.~Colours, data = slc_dford[ , c(comps_have, "Colours")], FUN = mean)

            colnames(centroidsOmniD)[2:ncol(centroidsOmniD)] <- paste0("meanPC", 1:length(comps_have))
            rownames(centroidsOmniD) <- centroidsOmniD$Colours
            #Get Euclidean distance between centroids
            cat(paste("Calculating OmniD Euclidean Distances in", slc, "\n"))
            centroiddistOmniD <- as.matrix(dist(centroidsOmniD[,2:ncol(centroidsOmniD)], method = "euclidean"))

            #Pump value into Slicedf
            Slicedf[slc , "Euclidean_Distance_P_NP_OmniD"] <- centroiddistOmniD["Progression", "No_Progression"]

        }


        aesthetic <- aes(x = PC1, y = PC2)
        p <- ggplot(slc_dford, aesthetic)
        p <- p + aes(col = Colours)

        if (!stratifybytreatment){
            groupcols <- setNames(as.character(c("#bcc2c2", "#00a6ff", "#c40000")), as.character(c("Censored", "No_Progression", "Progression")))
            groupcols_propplot <- groupcols
        } else {
            groupcols <- setNames(as.character(c("#bcc2c2", "#00a6ff", "#c40000")), as.character(c("Censored", "No_Progression", "Progression")))
            groupcols_propplot <- setNames(as.character(c("#bcc2c2", "#00a6ff", "#0a00cc", "#c40000")), as.character(c("Censored", "Under_Treatment_No_Progression", "Post_Treatment_No_Progression", "Progression")))
        }

        p <- p + scale_color_manual(values = groupcols)
        signifmeas <- "PERMANOVA"
        slc_pcatit <- paste(c(pcatit, (paste(signifmeas, "p <", Slicedf[slc , "Permanova"]))), collapse = "\n")
        distmessage <- paste0("Dissimilarity index = ", distmethod)
        slc_pcatit <- paste(c(slc_pcatit, distmessage), collapse = "\n")

        #p <- p + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))
        p <- p + theme(legend.position = "bottom")
        p <- p + geom_point(size = dotsize) + labs(x = xl, y = yl)
        if (!(is.null(forceaspectratio))){
            p <- p + theme(aspect.ratio = (1 / forceaspectratio))
        }

        if (TRUE){
            p <- p + geom_segment(aes(x = meanPC1, y = meanPC2, xend = PC1, yend = PC2, colour = Colours), data = centroiddf)
            #p <- p + geom_point(aes(x = meanPC1, y = meanPC2, colour = Colours), data = centroids, size = dotsize)
            p <- p + geom_point(aes(x = meanPC1, y = meanPC2), colour = "black", data = centroids, size = (dotsize * 3)) + geom_point(aes(x = meanPC1, y = meanPC2, colour = Colours), data = centroids, size = dotsize)
        }

        #p <- p + ggtitle(slc_pcatit)
        p <- p + labs(colour = slc)
        #p <- p + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
        if (grid == FALSE){
            p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }

        ##
        #Make secondary Ordination plot
        p2 <- ggplot(slc_dford, aesthetic)
        p2 <- p2 + aes(col = FollowUp_Time_days)

        p2 <- p2 + aes(shape = Colours)
        p2 <- add_shape_to_plot_safely(p = p2, shapevec = slc_dford$Colours, shapeby = "Colours", cdict = NULL)

        p2 <- p2 + scale_color_gradient(low = "blue", high = "red")

        #p2 <- p2 + theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size = 10))
        p2 <- p2 + theme(legend.position = "bottom")
        p2 <- p2 + geom_point(size = 2) + labs(x = xl, y = yl)
        if (!(is.null(forceaspectratio))){
            p2 <- p2 + theme(aspect.ratio = (1 / forceaspectratio))
        }

        #p2 <- p2 + ggtitle(slc_pcatit)

        p2 <- p2 + labs(shape = slc)
        #p2 <- p2 + guides(fill = guide_legend(nrow = 2, byrow = TRUE))


        if (grid == FALSE){
            p2 <- p2 + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        }

        finalordplot <- ggarrange(p, p2, ncol = 2, nrow = 1, common.legend = FALSE)
        finalordplot <- finalordplot + guides(fill = guide_legend(nrow = 3, byrow = TRUE))
        finalordplot <- annotate_figure(finalordplot, top = text_grob(slc_pcatit, color = "black", size = 8))

        SliceOrdinations[[slc]] <- finalordplot

    }

    #############################
    ## Get OR stats if required

    if (add_OR_asterisk){
        #Figure out which column to use
        if ("N_A" %in% currpt[ , "ORR"]){
            ORcolm <- "Study_Clin_Response"
        } else {
            ORcolm <- "ORR"
        }
        flog.info(paste("Using column", ORcolm, "as objective response data for adding asterisk."))

        dford_ORR <- dford
        colnames(dford_ORR)[which(colnames(dford_ORR) == ORcolm)] <- "ORRcategory"
        dford_ORR <- dford_ORR[ , c(comps_have, "ORRcategory")]

        ORR_stats_df <- data.frame(Measure = c("Permanova", "SumsOfSqs", "MeanSqs", "FModel", "R2", "Median_dist", "Mean_dist", "ANOSIM", "Reciprocal_Permanova", "Reciprocal_ANOSIM", "Euclidean_Distance_P_NP_2D", "Euclidean_Distance_P_NP_OmniD"), Value = NA)
        rownames(ORR_stats_df) <- ORR_stats_df$Measure

        set.seed(4140)
        #Get countmat with uncensored samples
        ValidSamples <- rownames(dford_ORR)[(dford_ORR[ , "ORRcategory"] != "N_A")]
        ORR_countmat <- countmat[ValidSamples, ]
        d <- vegdist(ORR_countmat, method = distmethod, na.rm = TRUE)

        ORR_pheno <- currpt[ValidSamples, ]
        colnames(ORR_pheno)[which(colnames(ORR_pheno) == ORcolm)] <- "ORRcategory"

        ORR_stats_df["Permanova", "Value"] <- vegan::adonis(as.formula(paste("d ~ ", "ORRcategory")), data = ORR_pheno, permutations = permanova_permutations)$aov.tab$`Pr(>F)`[1]

        ORR_stats_df["SumsOfSqs", "Value"] <- vegan::adonis(as.formula(paste("d ~ ", "ORRcategory")), data = ORR_pheno, permutations = permanova_permutations)$aov.tab$SumsOfSqs[1]

        ORR_stats_df["MeanSqs", "Value"] <- vegan::adonis(as.formula(paste("d ~ ", "ORRcategory")), data = ORR_pheno, permutations = permanova_permutations)$aov.tab$MeanSqs[1]

        ORR_stats_df["FModel", "Value"] <- vegan::adonis(as.formula(paste("d ~ ", "ORRcategory")), data = ORR_pheno, permutations = permanova_permutations)$aov.tab$F.Model[1]

        ORR_stats_df["R2", "Value"] <- vegan::adonis(as.formula(paste("d ~ ", "ORRcategory")), data = ORR_pheno, permutations = permanova_permutations)$aov.tab$R2[1]

        ORR_stats_df["ANOSIM", "Value"] <- anosim(ORR_countmat, as.factor(ORR_pheno[ , "ORRcategory"]), distance = distmethod, permutations = permanova_permutations)$signif

        #Get straight median dissimilarities, without reducing ordination
        curr_responders <- rownames(ORR_pheno)[which(ORR_pheno[ , "ORRcategory"] == "Responder")]
        curr_non_responders <- rownames(ORR_pheno)[which(ORR_pheno[ , "ORRcategory"] == "Non_Responder")]

        ORR_stats_df["Median_dist", "Value"] <- median(as.matrix(d)[curr_progressors, curr_no_progressors])
        ORR_stats_df["Mean_dist", "Value"] <- mean(as.matrix(d)[curr_progressors, curr_no_progressors])

        ORR_stats_df["Reciprocal_Permanova", "Value"] <- round((1 / ORR_stats_df["Permanova", "Value"]), 3)
        ORR_stats_df["Reciprocal_ANOSIM", "Value"] <- round((1 / ORR_stats_df["ANOSIM", "Value"]), 3)


        #Get Euclidean distance between centroids 2D
        centroids_ORR_2D <- aggregate(cbind(PC1, PC2) ~ ORRcategory, dford_ORR, mean)
        colnames(centroids_ORR_2D)[c(2, 3)] <- c("meanPC1", "meanPC2")
        rownames(centroids_ORR_2D) <- centroids_ORR_2D$ORRcategory

        centroiddist_ORR_2D <- as.matrix(dist(centroids_ORR_2D[,2:ncol(centroids_ORR_2D)], method = "euclidean"))

        #Pump value into Slicedf
        ORR_stats_df["Euclidean_Distance_P_NP_2D", "Value"] <- centroiddist_ORR_2D["Non_Responder", "Responder"]

        if (all(paste0("PC", 1:3) %in% colnames(dford_ORR))){
            centroidsOmniD <- aggregate(.~ORRcategory, data = dford_ORR[ , c(comps_have, "ORRcategory")], FUN = mean)
            colnames(centroidsOmniD)[2:ncol(centroidsOmniD)] <- paste0("meanPC", 1:length(comps_have))
            rownames(centroidsOmniD) <- centroidsOmniD$ORRcategory

            #Get Euclidean distance between centroids
            centroiddistOmniD <- as.matrix(dist(centroidsOmniD[,2:ncol(centroidsOmniD)], method = "euclidean"))

            #Pump value into Slicedf
            ORR_stats_df["Euclidean_Distance_P_NP_OmniD", "Value"] <- centroiddistOmniD["Non_Responder", "Responder"]
        }
    } else {
        ORR_stats_df <- NULL
    }

    ###########
    ## plot plots

    #Define useful plotting function
    plot_signifplot <- function(Slicedf = NULL, convert_days_to_months = TRUE, wantedcolname = NULL, loess_span = NULL, plottit = NULL, y_label = NULL, ORR_stats_df = NULL){

        curr_Slicedf <- Slicedf
        colnames(curr_Slicedf)[which(colnames(curr_Slicedf) == wantedcolname)] <- "Distance"

        if (convert_days_to_months){
            signif_plot <- ggplot(curr_Slicedf, aes(x = Months, y = Distance)) + geom_smooth(method = "loess", se = FALSE, span = loess_span) + geom_point(size = 1)
            xlab_name <- "Follow Up Months"
            signif_plot <- signif_plot + xlab(xlab_name) + ylab(y_label)
            xbreaks <- unique(c((curr_Slicedf$Months[seq(1, length(curr_Slicedf$Months), 3)]), max(curr_Slicedf$Months)))
            xlabels <- as.character(xbreaks)
            ORRtime <- 90 / 30
        } else {
            signif_plot <- ggplot(curr_Slicedf, aes(x = Days, y = Distance)) + geom_smooth(method = "loess", se = FALSE, span = loess_span) + geom_point(size = 1)
            xlab_name <- "Follow Up Days"
            signif_plot <- signif_plot + xlab(xlab_name) + ylab(y_label)
            xbreaks <- unique(c((curr_Slicedf$Days[seq(1, length(curr_Slicedf$Days), 3)]), max(curr_Slicedf$Days)))
            xlabels <- as.character(xbreaks)
            ORRtime <- 90
        }

        if (!is.null(ORR_stats_df)){
            wantedval <- ORR_stats_df[wantedcolname , "Value"]
            signif_plot <- signif_plot + geom_point(aes(x = ORRtime, y = wantedval), shape = 8, colour = "blue")
        }

        signif_plot <- signif_plot + scale_x_continuous(breaks = xbreaks, labels = xlabels)

        signif_plot <- signif_plot + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
        #signif_plot <- signif_plot + geom_hline(yintercept = (1/0.05), linetype = "dashed", color = "red")

        signif_plot <- signif_plot + theme(title = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"), axis.title.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"), axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"))

        signif_plot <- signif_plot + ggtitle(plottit)

        return(signif_plot)

    }

    if (!is.null(crop_to_time)){
        if (convert_days_to_months){
            Slicedf <- subset(Slicedf, Months <= crop_to_time)
            #If there is anything smaller than 3 months, like 2 or 1, cut that.
            Slicedf <- subset(Slicedf, Months >= 3)
        } else {
            Slicedf <- subset(Slicedf, Days <= crop_to_time)
        }
    }

    #Plot permanova significance
    if (convert_days_to_months){
        signif_plot <- ggplot(Slicedf, aes(x = Months, y = Reciprocal_Permanova)) + geom_smooth(method = "loess", se = FALSE, span = loess_span) + geom_point(size = 1)
        xlab_name <- "Follow Up Months"
        signif_plot <- signif_plot + xlab(xlab_name) + ylab(paste("1", "/", signifmeas))
        xbreaks <- unique(c((Slicedf$Months[seq(1, length(Slicedf$Months), 3)]), max(Slicedf$Months)))
        xlabels <- as.character(xbreaks)
        ORRtime <- 90 / 30
    } else {
        signif_plot <- ggplot(Slicedf, aes(x = Days, y = Reciprocal_Permanova)) + geom_smooth(method = "loess", se = FALSE, span = loess_span) + geom_point(size = 1)
        xlab_name <- "Follow Up Days"
        signif_plot <- signif_plot + xlab(xlab_name) + ylab(paste("1", "/", signifmeas))
        xbreaks <- unique(c((Slicedf$Days[seq(1, length(Slicedf$Days), 3)]), max(Slicedf$Days)))
        xlabels <- as.character(xbreaks)
        ORRtime <- 90
    }

    signif_plot <- signif_plot + scale_x_continuous(breaks = xbreaks, labels = xlabels)

    signif_plot <- signif_plot + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
    signif_plot <- signif_plot + geom_hline(yintercept = (1/0.05), linetype = "dashed", color = "black")

    if (max(Slicedf$Reciprocal_Permanova) > (1 / 0.01)){
        signif_plot <- signif_plot + geom_hline(yintercept = (1 / 0.01), linetype = "dashed", color = "black")
    }


    if (!is.null(ORR_stats_df)){
        wantedval <- ORR_stats_df["Reciprocal_Permanova" , "Value"]
        signif_plot <- signif_plot + geom_point(aes(x = ORRtime, y = wantedval), shape = 8, colour = "blue")
    }

    dist_plot_2D <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Euclidean_Distance_P_NP_2D", loess_span = loess_span, plottit = "Distance between P/NP Centroids in 2D", y_label = "Distance between P/NP Centroids in 2D")

    if ("Euclidean_Distance_P_NP_OmniD" %in% colnames(Slicedf)){
        dist_plot_OmniD <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Euclidean_Distance_P_NP_OmniD", loess_span = loess_span, plottit = paste("Distance between P/NP Centroids in", length(comps_have), "Principal Components"), y_label = "Distance between P/NP Centroids in OmniD", ORR_stats_df = ORR_stats_df)
    }


    FModel_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "FModel", loess_span = loess_span, plottit = "Permanova F-Model plot", y_label = "Permanova F-Model")

    SumsOfSqs_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "SumsOfSqs", loess_span = loess_span, plottit = "Permanova SumsOfSqs plot", y_label = "Permanova SumsOfSqs")

    MeanSqs_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "MeanSqs", loess_span = loess_span, plottit = "Permanova MeanSqs plot", y_label = "Permanova MeanSqs")

    R2_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "R2", loess_span = loess_span, plottit = "Permanova R2 plot", y_label = "Permanova R2")

    ANOSIM_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Reciprocal_ANOSIM", loess_span = loess_span, plottit = "ANOSIM plot", y_label = "1 / ANOSIM p-value")

    Median_dist_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Median_dist", loess_span = loess_span, plottit = "Median_dist plot", y_label = "Median_dist")

    Mean_dist_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Mean_dist", loess_span = loess_span, plottit = "Mean_dist plot", y_label = "Mean_dist")

    #Min_dist_plot <- plot_signifplot(Slicedf = Slicedf, convert_days_to_months = convert_days_to_months, wantedcolname = "Min_dist", loess_span = loess_span, plottit = "Min_dist plot", y_label = "Min_dist")

    #Plot the proportions

    if (stratifybytreatment){
        longdat <- Slicedf %>% gather(Patient_Type, Number_of_Samples, c(Censored, Under_Treatment_No_Progression, Post_Treatment_No_Progression, Progression))
        longdat$Patient_Type <- factor(longdat$Patient_Type, levels = rev(c("Progression", "Under_Treatment_No_Progression", "Post_Treatment_No_Progression", "Censored")))

    } else {
        longdat <- Slicedf %>% gather(Patient_Type, Number_of_Samples, c(Censored, No_Progression, Progression))
        longdat$Patient_Type <- factor(longdat$Patient_Type, levels = rev(c("Progression", "No_Progression", "Censored")))
    }

    if (convert_days_to_months){
        proportionsplot <- ggplot(longdat, aes(x = Months, y = Number_of_Samples, fill = Patient_Type)) + geom_area()
    } else {
        proportionsplot <- ggplot(longdat, aes(x = Days, y = Number_of_Samples, fill = Patient_Type)) + geom_area()
    }

    #groupcols <- setNames(as.character(c("#bcc2c2", "#c40000", "#00a6ff")), as.character(c("Censored", "No_Progression", "Progression")))
    proportionsplot <- proportionsplot + scale_fill_manual(values = groupcols_propplot)
    proportionsplot <- proportionsplot + scale_x_continuous(breaks = xbreaks, labels = xlabels)
    proportionsplot <- proportionsplot + xlab(xlab_name) + ylab("Number of Samples")

    proportionsplot <- proportionsplot + theme(legend.position =  "bottom", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))

    if ("Patients_under_treatment" %in% colnames(Slicedf)){
        if (convert_days_to_months){
            proportionsplot <- proportionsplot + geom_line(aes(x = Months, y = Patients_under_treatment))
        } else {
            proportionsplot <- proportionsplot + geom_line(aes(x = Days, y = Patients_under_treatment))
        }
    }

    if (FALSE){
        wantedplots <- c("dist_plot_2D", "dist_plot_10D", "dist_plot_60D", "FModel_plot", "SumsOfSqs_plot", "MeanSqs_plot", "R2_plot", "ANOSIM_plot", "Median_dist_plot", "Mean_dist_plot")
    } else {
        wantedplots <- c("dist_plot_2D", "dist_plot_OmniD", "FModel_plot", "SumsOfSqs_plot", "MeanSqs_plot", "R2_plot", "ANOSIM_plot", "Median_dist_plot", "Mean_dist_plot")
    }

    gvec <- lapply(wantedplots, function(x){ ggarrange(proportionsplot, get(x), signif_plot, ncol = 1, nrow = 3, heights = c(1, 1, 1), vjust = 1, hjust = -1) })

    names(gvec) <- wantedplots

    if (algorithm != "tUMAP"){
        gvec[[length(gvec) + 1]] <- varplot
        names(gvec)[length(gvec)] <- "Variance_by_Principal_Component_plot"
    }

    if (return_Ordinations){
        gvec <- c(gvec, SliceOrdinations)
    }

    if (return_metadata_slices){
        gvec$metadata_slices <- currpt
    }

    if (return_Slice_df){
        gvec$Slicedf <- Slicedf
    }

    return(gvec)

}
