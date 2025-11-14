# curtailment 1.0.0

 - Now using `plot` for singlearmDesign and find2stagedesigns (remaining designs to follow, then drawDiagram may become deprecated)
 - "Breaking" updates to argument names nmin, nmax -- now using argument n.max, vector of length 2.

# curtailment 0.3.0

find2stageDesigns now orders by increasing Nmax, row numbers removed.

# curtailment 0.2.7

When plotting single-arm designs and general plots (drawDiagramGeneric), added argument benefit.stop to allow diagrams to show no benefit stopping 

# curtailment 0.2.6

Changed output of find2stageDesigns

# curtailment 0.2.5

* Tidied output of find2stageDesigns

# curtailment 0.2.4

* Fixed bug in simonEfficacy

# curtailment 0.2.3

* Sped up Mander & Thompson design search (command find2stageDesigns)

# curtailment 0.2.2

* Big fix: minstop argument in singlearmDesign

# curtailment 0.2.1

* Bug fixes

# curtailment 0.2.0

* Added minstop argument across multiple functions, specifying the minimum sample size for which early stopping is permitted.

* Updated the find2stageDesigns function to allow a maximum permitted conditional power.

# curtailment 0.1.1

* Changed the description, which now does not start with the package name.

* Reduced the length of the title to less than 65 characters.

* Added \value (including output class) to the following .Rd files:
  * drawDiagramGeneric.Rd
  * findNSCdesigns.Rd
  * findSCdesigns.Rd
  * plotByWeight.Rd

# curtailment 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added a `README.md` file to describe the package
* Updated all functions to comply with R CMD check.
