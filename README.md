# OceanOmics 'super tree'

For visualisation purposes we need a large tree that contains all Australian marine mammals. For subsequent publications, presentations, 
visualisations etc. we can then always 'just' subset that large tree.

This repository contains code that pulls in several public trees, 
merges them, 
and then treats that large tree as a backbone to insert, at the genus-level, species we do not have in the phylogeny but are Australian marine vertebrates.


# Caveats

- The branch lengths in the super tree are meaningless as the contributing trees do things differently.

- The Australian marine vertebrates are 'randomly' placed into their genus (see [here](https://cran.r-project.org/web/packages/RRphylo/vignettes/Tree-Manipulation.html) for the method),
so their placement is a bit random, and obviously their branch length is extra-meaningless.

# Files

marine_verts.txt - an OceanOmics-internal curated list of Australian marine vertebrate species. Probably has some errors.

big_tree.Actinopterygii.Chondrichtyes.tre - the Betancur-R tree and the Stein tree merged, nothing else.

big_tree.Actinopterygii.Chondrichtyes.withExtraSpecies.tre - the above tree with as many Australian marine vertebrates added as possible. Where Genus was not present in the tree, the species was not added.

-- more to come here

.Renviron - this just changes the download SSL package to use openSSL to work with our internal VPN.

# To subset

```
# load the tree
t_big <- ape::read.tree("big_tree.Actinopterygii.Chondrichtyes.withMarineVerts.tre")

# make sure your list of species has _ in their names instead of spaces
your_list_of_species <- gsub(' ', '_', all_species)

# check which of those species are in the tree
species_in_tree <- intersect(species_to_keep, t_big$tip.label)

# and subset the tree
t2 <- ape::keep.tip(t_big, species_in_tree)
```

# Data Sources/Citations

Actinopterygii: Betancur-R, R., Wiley, E.O., Arratia, G. et al. Phylogenetic classification of bony fishes. BMC Evol Biol 17, 162 (2017). https://doi.org/10.1186/s12862-017-0958-3 Additional File 2

Actinopterygii, second tree: via the fishtreeoflife R package, Chang J, Rabosky DL, Smith SA, Alfaro ME (2019). “An R package and online resource for macroevolutionary studies using the ray-finned fish tree of life.” Methods in Ecology and Evolution, 10(7), 1118-1124. ISSN 2041-210X, doi:10.1111/2041-210x.13182.

Rabosky DL, Chang J, Title PO, Cowman PF, Sallan L, Friedman M, Kaschner K, Garilao C, Near TJ, Coll M, Alfaro ME (2018). “An inverse latitudinal gradient in speciation rate for marine fishes.” Nature, 559(7714), 392-395. doi:10.1038/s41586-018-0273-1.

Chondrichtyes: Global priorities for conserving the evolutionary history of sharks, rays, and chimaeras R.W. Stein, C.G. Mull, T.S. Kuhn, N.C. Aschliman, L.N.K. Davidson, J.B. Joy, G.J. Smith, N.K. Dulvy, A.O. Mooers Nature Ecolgy and Evolution, 2. 
from https://vertlife.org/sharktree/downloads/ 610.tree.10Cal.RAxML.BS.nex
 
Aves:

(Marine) Mammalia:

Seasnakes:
