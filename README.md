# GCM_plots
This Repository contains python functions and notebooks designed to replicated the original IDL code for generating outputs from RM-GCM


-> Each notebook ~replicates an old IDL script, with an associated .py function

-> Each notebook (hopefully) contains detailed instructions for what each necessary input is

-> To use: first run the FORTRAN FILE FINISHER
    this reads and converts relevant fort.* files and saves into a pandas data frame
-> the E_IGCM notebook can make lat v. lon plots of all the things! (winds,temps, outgoing radiations) 

Note: since git sees a "change" to notebooks when you just run them and don't actually change them, you can use git stash to "hide" that from git before you do a git pull
