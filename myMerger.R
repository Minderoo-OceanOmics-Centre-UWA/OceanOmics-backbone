
my_merger <- function (backbone = NULL, data, source.tree = NULL, age.offset = NULL, 
                       tip.ages = NULL, node.ages = NULL, plot = TRUE, filename = NULL) 
{
  misspacks <- sapply(c("manipulate", "scales"), requireNamespace, 
                      quietly = TRUE)
  if (any(!misspacks)) {
    stop("The following package/s are needed for this function to work, please install it/them:\n ", 
         paste(names(misspacks)[which(!misspacks)], collapse = ", "), 
         call. = FALSE)
  }
  dat <- data
  if (is.null(backbone)) {
    if (dat[1, ]$reference %in% dat$bind) 
      stop("Please indicate as first row of data the first pair of tips to build the tree on.")
    startsp <- unlist(dat[1, 1:2])
    dat <- dat[-1, ]
    if (any(dat$reference == startsp[2])) {
      dat$poly[which(dat$reference == startsp[2])] <- TRUE
      dat$reference[which(dat$reference == startsp[2])] <- paste(startsp, 
                                                                 collapse = "-")
    }
    tree <- list()
    tree$edge <- matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE)
    tree$tip.label <- startsp
    tree$Nnode <- 1
    tree$edge.length <- rep(1, 2)
    class(tree) <- "phylo"
  }
  else tree <- backbone
  tree2 <- source.tree
  H <- max(diag(vcv(tree)))
  ages <- H - diag(vcv(tree))
  if (!is.null(age.offset) && age.offset < 0) {
    ages <- ages + abs(age.offset)
    Hset <- H + abs(age.offset)
  }
  else Hset <- H
  if (!all(colnames(dat) %in% c("bind", "reference", "poly"))) {
    if (any(is.na(as.logical(dat[, 3])))) 
      stop("Check columns order: it should be 'bind', 'reference', 'poly'")
    warning("Colnames not matching: columns assumed to be ordered as 'bind','reference','poly\n'", 
            immediate. = TRUE)
    colnames(dat) <- c("bind", "reference", "poly")
  }
  if (!is.logical(dat$poly)) 
    dat$poly <- as.logical(dat$poly)
  if (any(apply(dat, 1, function(k) (!grepl("-", k[2])) & !grepl("Genus", 
                                                                 k[2]) & !grepl("Clade", k[2]) & as.logical(k[3])))) {
    warning("Found poly=TRUE at binding two tips, automatically changed to FALSE\n", 
            immediate. = TRUE)
    lapply(1:nrow(dat), function(k) {
      if ((!grepl("-", dat[k, 2])) & !grepl("Genus", k[k, 
                                                       2]) & !grepl("Clade", k[k, 2]) & dat[k, 3]) {
        dat[k, 3] <<- FALSE
        message(paste("  Please check", dat[k, 1], "&", 
                      dat[k, 2], "in your dataset", sep = " "))
      }
    })
  }
  if (any(grepl("Genus", dat$bind))) {
    if (is.null(tree2)) 
      stop("Please provide the source.tree")
    sapply(1:nrow(dat), function(k) {
      if (grepl("Genus", dat[k, ]$bind)) {
        genref <- trimws(gsub("Genus ", "", dat[k, ]$bind))
        getgenref <- getGenus(tree2, genref)
        if (getgenref[2] == 1) 
          dat[k, ]$bind <<- grep(genref, tree2$tip.label, 
                                 value = TRUE)
        else dat[k, ]$bind <<- paste(tips(tree2, getgenref[, 
                                                           3])[c(1, length(tips(tree2, getgenref[, 3])))], 
                                     collapse = "-")
      }
    })
  }
  if (any(grepl("Clade", dat$bind))) {
    sapply(1:nrow(dat), function(k) {
      if (grepl("Clade", dat[k, ]$bind)) {
        if (!gsub("Clade ", "", dat[k, ]$bind) %in% tree2$node.label) 
          stop("Required node.label not indicated on the source.tree")
        claref <- tips(tree2, Ntip(tree2) + which(tree2$node.label == 
                                                    gsub("Clade ", "", dat[k, ]$bind)))
        dat[k, ]$bind <<- paste(claref[c(1, length(claref))], 
                                collapse = "-")
      }
    })
  }
  dat <- data.frame(dat, bind.type = sapply(strsplit(dat[, 
                                                         1], "-"), length))
  if (any(dat$bind.type == 2) & is.null(tree2)) 
    stop("Please provide the source.tree")
  bind.all <- unlist(sapply(1:nrow(dat), function(x) {
    ifelse(dat[x, 4] == 2, des <- tips(tree2, getMRCA(tree2, 
                                                      strsplit(dat[x, 1], "-")[[1]])), des <- dat[x, 1])
    des
  }))
  if (any(duplicated(bind.all))) 
    stop(paste(paste(bind.all[duplicated(bind.all)], collapse = ", "), 
               "names duplicated in supplied tips"))
  if (all(dat$bind.type == 1) & (!is.null(tree2))) 
    tree2 <- NULL
  dat$bind.tips <- dat$bind
  if (any(dat$bind.type == 2)) 
    dat[which(dat$bind.type == 2), ]$bind.tips <- lapply(dat$bind[which(dat$bind.type == 
                                                                          2)], function(x) tips(tree2, getMRCA(tree2, strsplit(x, 
                                                                                                                               "-")[[1]])))
  if (!all(unlist(dat[which(dat$bind.type == 2), ]$bind.tips) %in% 
           tree2$tip.label)) {
    stop(paste(paste(unlist(dat[which(dat$bind.type == 2), 
    ]$bind.tips)[which(unlist(dat[which(dat$bind.type == 
                                          2), ]$bind.tips) %in% tree2$tip.label)], collapse = ", "), 
    "not in source.tree"))
  }
  if (any(dat$bind %in% tree$tip.label)) {
    warning(paste(paste(dat[which(dat$bind %in% tree$tip.label), 
                            1], collapse = ", "), "removed from the backbone tree\n"), 
            immediate. = TRUE)
    tree <- drop.tip(tree, dat[which(dat$bind %in% tree$tip.label), 
                               1])
  }
  if (any(grepl("Genus", dat$reference))) {
    dat$bindgen <- sapply(dat$bind.tips, function(j) unique(sapply(j, 
                                                                   function(x) strsplit(x, "_")[[1]][1])))
    sapply(1:nrow(dat), function(k) {
      if (grepl("Genus", dat[k, ]$reference)) {
        genref <- trimws(gsub("Genus ", "", dat[k, ]$reference))
        getgenref <- try(getGenus(tree, genref), silent = TRUE)
        if (!inherits(getgenref, "try-error")) {
          if (getgenref[2] == 1) 
            gentree <- grep(genref, tree$tip.label, value = TRUE)
          else gentree <- paste(tips(tree, getgenref[, 
                                                     3])[c(1, length(tips(tree, getgenref[, 3])))], 
                                collapse = "-")
          genbind <- dat$bindgen[[k]]
          dat[k, ]$reference <<- gentree
          if (genref %in% genbind) {
            if (length(genbind) > 1 & any(sapply(dat$bindgen, 
                                                 function(w) all(w %in% genref)))) 
              dat[k, ]$reference <<- paste(c(gentree, 
                                             sapply(dat[which(sapply(dat$bindgen, 
                                                                     function(w) all(w %in% genref))), ]$bind.tips, 
                                                    "[[", 1)), collapse = "-")
          }
          else {
            if (any(sapply(dat$bindgen, function(w) any(w %in% 
                                                        genref)))) 
              dat[k, ]$reference <<- paste(c(gentree, 
                                             sapply(dat[which(sapply(dat$bindgen, 
                                                                     function(w) any(w %in% genref))), ]$bind.tips, 
                                                    "[[", 1)), collapse = "-")
          }
        }
        else stop(paste("Genus", genref, "missing from the backbone tree"))
      }
    })
    dat$bindgen <- NULL
  }
  if (any(grepl("Clade", dat$reference))) {
    sapply(1:nrow(dat), function(k) {
      if (grepl("Clade", dat[k, ]$reference)) {
        if (!gsub("Clade ", "", dat[k, ]$reference) %in% 
            tree$node.label) 
          stop("Required node.label not indicated on the tree")
        claref <- tips(tree, Ntip(tree) + which(tree$node.label == 
                                                  gsub("Clade ", "", dat[k, ]$reference)))
        dat[k, ]$reference <<- paste(claref[c(1, length(claref))], 
                                     collapse = "-")
      }
    })
  }
  tab.ref <- table(dat$reference)
  if (any(tab.ref > 1)) {
    for (j in 1:length(which(tab.ref > 1))) {
      ref.mult <- dat[which(dat$reference == names(which(tab.ref > 
                                                           1)[j])), ]
      if (any(ref.mult$poly)) 
        ref.mult <- ref.mult[which(!ref.mult$poly), ]
      ref.mult$reference[-1] <- paste(strsplit(ref.mult$reference[1], 
                                               "-")[[1]][1], strsplit(ref.mult$bind[1], "-")[[1]][1], 
                                      sep = "-")
      ref.mult$poly[-1] <- TRUE
      dat[match(ref.mult$bind, dat$bind), ] <- ref.mult
    }
  }
  if (!is.null(tree2)) {
    dat$MRCAbind <- NA
    dat$MRCAbind[which(dat$bind.type == 2)] <- sapply(dat[which(dat$bind.type == 
                                                                  2), ]$bind.tips, function(x) getMRCA(tree2, x))
    if (any(dat$bind.type == 2)) {
      remt <- unlist(sapply(dat[which(dat$bind.type == 
                                        2), ]$MRCAbind, function(x) tree$tip.label[which(tree$tip.label %in% 
                                                                                           tips(tree2, x))]))
      if (length(remt) > 0) {
        tree <- drop.tip(tree, remt)
        ages <- ages[-match(remt, names(ages))]
        warning(paste(paste(remt, collapse = ", "), "already on the source tree: removed from the backbone tree\n"), 
                immediate. = TRUE)
      }
    }
  }
  refs <- strsplit(dat$reference, "-")
  dat$ref.tree1 <- NA
  ref.tree <- lapply(refs, function(x) {
    if (all(x %in% tree$tip.label)) 
      NA
    else x[which(!x %in% tree$tip.label)]
  })
  dat[which(is.na(ref.tree)), ][order(dat[which(is.na(ref.tree)), 
  ]$bind.type), ]$ref.tree1 <- seq(1, length(which(is.na(ref.tree))))
  if (any(!is.na(ref.tree))) {
    dat.new <- dat[which(!is.na(ref.tree)), ]
    dat.new$ref.tree1 <- ref.tree[-which(is.na(ref.tree))]
    if (any(unlist(dat.new$ref.tree1) %in% unlist(dat.new$bind.tips))) {
      while (nrow(dat.new) > 0) {
        outs <- which(!sapply(dat.new$ref.tree1, function(w) any(w %in% 
                                                                   unlist(dat.new$bind.tips))))
        if (length(outs) < 1) 
          stop("Recursive species attachment: check rows ", 
               paste(rownames(dat.new), collapse = ", "), 
               " in data")
        dat[match(dat.new[outs, 1], dat[, 1]), ]$ref.tree1 <- max(dat$ref.tree1, 
                                                                  na.rm = TRUE) + 1:length(outs)
        dat.new <- dat.new[-outs, ]
      }
    }
    else {
      outs <- which(!sapply(dat.new$ref.tree1, function(w) any(w %in% 
                                                                 unlist(dat.new$bind.tips))))
      dat[match(dat.new[outs, 1], dat[, 1]), ]$ref.tree1 <- max(dat$ref.tree1, 
                                                                na.rm = TRUE) + 1:length(which(!dat.new$ref.tree1 %in% 
                                                                                                 dat.new[, 1]))
    }
  }
  dat <- dat[order(dat$ref.tree1), ]
  pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
  for (k in 1:nrow(dat)) {
    setTxtProgressBar(pb, k)
    if (length(strsplit(dat$reference, "-")[[k]]) > 1) 
      where.ref <- getMRCA(tree, strsplit(dat$reference, 
                                          "-")[[k]])
    else where.ref <- which(tree$tip.label == strsplit(dat$reference, 
                                                       "-")[[k]])
    print(strsplit(dat$reference, 
                   "-")[[k]])
    print(where.ref)
    print(Ntip(tree))
    if (where.ref != (Ntip(tree) + 1)) 
      br.len <- tree$edge.length[which(tree$edge[, 2] == 
                                         where.ref)]
    else {
      if (is.null(tree$root.edge) || tree$root.edge == 
          0) 
        tree$root.edge <- mean(tree$edge.length)
      br.len <- tree$root.edge
    }
    print(dat$bind.type[k])
    if (dat$bind.type[k] == 1) {
      if (dat$poly[k]) {
        pos.ref <- 0
        br.len <- max(tree$edge.length[which(tree$edge[, 
                                                       1] %in% where.ref)])
      }
      else pos.ref <- br.len/2

      tree <- bind.tip(tree, dat$bind[k], where.ref, position = pos.ref, 
                       edge.length = br.len/2)
    }
    else {
      cla <- extract.clade(tree2, dat$MRCAbind[k])
      if (dat$poly[k]) {
        pos.ref <- 0
        if (where.ref == (Ntip(tree) + 1) && (max(diag(vcv(cla))) + 
                                              max(diag(vcv(cla)))/10) > H) 
          cla <- rescaleRR(cla, height = (H + max(diag(vcv(cla)))/10))
        else if (where.ref != (Ntip(tree) + 1) & (max(diag(vcv(cla))) + 
                                                  max(diag(vcv(cla)))/10) > (H - nodeHeights(tree)[which(tree$edge[, 
                                                                                                                   2] == where.ref), 2])) 
          cla <- rescaleRR(cla, height = (H - nodeHeights(tree)[which(tree$edge[, 
                                                                                2] == where.ref), 2] + max(diag(vcv(cla)))/10))
        cla$root.edge <- max(diag(vcv(cla)))/10
      }
      else {
        if (where.ref == (Ntip(tree) + 1)) 
          pos.ref <- br.len/2 + (max(diag(vcv(cla))) - 
                                   max(diag(vcv(tree))))
        else {
          Hcla <- max(diag(vcv(cla))) + br.len/2
          if ((H - nodeHeights(tree)[which(tree$edge[, 
                                                     2] == where.ref), 2]) < H/1000) {
            cla <- rescaleRR(cla, height = br.len/4)
            pos.ref <- br.len - br.len/4
          }
          else {
            if (Hcla > (H - nodeHeights(tree)[which(tree$edge[, 
                                                              2] == where.ref), 1])) {
              cla <- rescaleRR(cla, height = (H - nodeHeights(tree)[which(tree$edge[, 
                                                                                    2] == where.ref), 2]))
              pos.ref <- br.len/2
            }
            else if (Hcla > (H - nodeHeights(tree)[which(tree$edge[, 
                                                                   2] == where.ref), 2])) 
              pos.ref <- (max(diag(vcv(cla))) + br.len/2) - 
                (H - nodeHeights(tree)[which(tree$edge[, 
                                                       2] == where.ref), 2])
            else pos.ref <- br.len/2
          }
        }
        cla$root.edge <- br.len/2
      }
      tree <- bind.tree(tree, cla, where = where.ref, position = pos.ref)
    }
  }
  close(pb)
  message("Binding done!")
  if (is.null(backbone)) {
    tree$edge.length <- rep(1, length(tree$edge.length))
    if (!is.null(tip.ages)) 
      tree <- rescaleRR(tree, height = tip.ages[which.max(tip.ages)] + 
                          tip.ages[which.max(tip.ages)]/10)
  }
  trycal <- try({
    if (is.null(tip.ages)) {
      tip.ages <- rep(0, length(tree$tip.label))
      names(tip.ages) <- tree$tip.label
      tip.ages[match(names(ages), names(tip.ages))] <- ages
    }
    else {
      if (any(!names(tip.ages) %in% tree$tip.label)) {
        warning(paste(paste(names(tip.ages)[which(!names(tip.ages) %in% 
                                                    tree$tip.label)], collapse = ", "), "not on the final tree: removed from the vector of tip.ages\n"), 
                immediate. = TRUE)
        tip.ages <- tip.ages[which(names(tip.ages) %in% 
                                     tree$tip.label)]
      }
      if (!all(names(ages) %in% names(tip.ages))) 
        tip.ages <- c(ages[which(!names(ages) %in% names(tip.ages))], 
                      tip.ages)
      if (!all(tree$tip.label %in% names(tip.ages))) {
        tip.add <- rep(0, length(which(!tree$tip.label %in% 
                                         names(tip.ages))))
        names(tip.add) <- tree$tip.label[which(!tree$tip.label %in% 
                                                 names(tip.ages))]
        tip.ages <- c(tip.ages, tip.add)
      }
    }
    if (!is.null(node.ages)) {
      names(node.ages) <- sapply(names(node.ages), function(x) {
        getMRCA(tree, unlist(strsplit(x, "-")))
      })
    }
    else node.ages <- c()
    if (!is.null(backbone) & !getMRCA(tree, names(ages)) %in% 
        names(node.ages)) 
      node.ages <- setNames(c(node.ages, Hset), c(names(node.ages), 
                                                  getMRCA(tree, names(ages))))
    if (any(dat$bind.type == 2)) {
      ages.fix <- unlist(lapply(dat$MRCAbind[which(dat$bind.type == 
                                                     2)], function(x) {
                                                       des <- c(x, getDescendants(tree2, x)[which(getDescendants(tree2, 
                                                                                                                 x) > Ntip(tree2))])
                                                       dn <- dist.nodes(tree2)
                                                       dndes <- max(diag(vcv(tree2))) - dn[which(rownames(dn) == 
                                                                                                   (Ntip(tree2) + 1)), match(des, rownames(dn))]
                                                       names(dndes) <- des
                                                       dndes
                                                     }))
      if (!is.null(age.offset) && age.offset > 0) 
        ages.fix <- ages.fix + age.offset
      names(ages.fix) <- sapply(1:length(ages.fix), function(x) {
        nodex <- getMRCA(tree, tips(tree2, as.numeric(names(ages.fix)[x])))
        tipx.ages <- tip.ages[match(tips(tree, nodex), 
                                    names(tip.ages))]
        if (any(tipx.ages > ages.fix[x])) 
          NA
        else nodex
      })
      ages.fix <- ages.fix[which(!is.na(names(ages.fix)))]
      if (any(!names(ages.fix) %in% names(node.ages))) 
        node.ages <- c(ages.fix[which(!names(ages.fix) %in% 
                                        names(node.ages))], node.ages)
    }
    if (max(diag(vcv(tree))) > H && (!(Ntip(tree) + 1) %in% 
                                     names(node.ages))) 
      warning(paste("Root age not indicated: the tree root arbitrarily set at\n", 
                    round(max(diag(vcv(tree))), 2)), immediate. = TRUE)
    tree.plot <- tree.final <- scaleTree(tree, node.ages = node.ages, 
                                         tip.ages = tip.ages)
    message("Age Calibration done!")
  }, silent = TRUE)
  if (inherits(trycal, "try-error")) {
    warning("Age Calibration failed with the error: \n", 
            trycal[1], "\n Returning the uncalibrated version of the tree", 
            call. = FALSE)
    tree.plot <- tree.final <- tree
  }
  if (plot) {
    if (any(dat$bind.type == 2)) 
      cla.plot <- lapply(dat$MRCAbind[which(dat$bind.type == 
                                              2)], function(x) c(getMRCA(tree.plot, tips(tree2, 
                                                                                         x)), getDescendants(tree.plot, getMRCA(tree.plot, 
                                                                                                                                tips(tree2, x)))))
    else cla.plot <- c()
    all.plot <- c(which(tree.plot$tip.label %in% dat$bind[which(dat$bind.type == 
                                                                  1)]), unlist(cla.plot))
    colo <- rep((scales::hue_pal())(2)[2], nrow(tree.plot$edge))
    colo[match(all.plot, tree.plot$edge[, 2])] <- (scales::hue_pal())(2)[1]
    names(colo) <- tree.plot$edge[, 2]
    colo.tips <- colo[which(as.numeric(names(colo)) <= Ntip(tree.plot))]
    tree.plot$tip.label[which(!tree.plot$tip.label %in% unique(c(unlist(strsplit(dat$bind, 
                                                                                 "-")), unlist(strsplit(dat$reference, "-")))))] <- " "
    names(colo.tips) <- tree.plot$tip.label[as.numeric(names(colo.tips))]
    if (!is.null(filename)) {
      pdf(file = paste(filename, ".pdf", sep = ""))
      if (Ntip(tree.plot) < 100) 
        plot(tree.plot, edge.color = colo, cex = 0.6, 
             tip.color = colo.tips[match(tree.plot$tip.label, 
                                         names(colo.tips))])
      else plot(tree.plot, edge.color = colo, type = "fan", 
                cex = 0.6, tip.color = colo[which(tree.plot$edge[, 
                                                                 2] <= Ntip(tree.plot))])
      dev.off()
    }
    clades.plot <- lapply(1:nrow(dat), function(x) {
      MRCAplot <- getMRCA(tree.final, c(unlist(dat$bind.tips[x]), 
                                        unlist(strsplit(dat$reference[x], "-"))))
      cla <- extract.clade(tree.final, getMRCA(tree.final, 
                                               c(unlist(dat$bind.tips[x]), unlist(strsplit(dat$reference[x], 
                                                                                           "-")))))
      colo.cla <- colo[which(names(colo) %in% getDescendants(tree.final, 
                                                             MRCAplot))]
      nam.colo <- as.numeric(names(colo.cla))
      nam.colo[which(nam.colo <= Ntip(tree.final))] <- which(cla$tip.label %in% 
                                                               tree.final$tip.label[nam.colo[which(nam.colo <= 
                                                                                                     Ntip(tree.final))]])
      if (any(nam.colo > Ntip(tree.final))) 
        nam.colo[which(nam.colo > Ntip(tree.final))] <- sapply(nam.colo[which(nam.colo > 
                                                                                Ntip(tree.final))], function(r) getMRCA(cla, 
                                                                                                                        tips(tree.final, r)))
      names(colo.cla) <- nam.colo
      colo.cla <- colo.cla[match(cla$edge[, 2], names(colo.cla))]
      list(cla, colo.cla)
    })
    names(clades.plot) <- dat[, 1]
    taxon <- manipulate::picker(as.list(dat[, 1]))
    manipulate::manipulate(plot(clades.plot[[which(names(clades.plot) == 
                                                     taxon)]][[1]], edge.color = clades.plot[[which(names(clades.plot) == 
                                                                                                      taxon)]][[2]], cex = 0.6), taxon = taxon)
  }
  return(tree.final)
}