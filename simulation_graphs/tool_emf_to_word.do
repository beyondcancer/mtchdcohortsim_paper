// Define the folder containing the .emf files

local folder "$resultsDir"

// List all the .emf files in the folder
local files : dir "`folder'" files "*ss.emf"

// Create a new Word document
putdocx begin

// Loop through each .emf file and insert it into the Word document
foreach file of local files {
    
    // Insert the .emf graph into the Word document
    putdocx paragraph
    putdocx image "`folder'/`file'"
    
    // Optionally, you can add titles or additional text for each graph
    putdocx text ("Graph: `file'")
}

// Save and close the Word document
putdocx save "$resultsDir/stata_graphs.docx", replace

putdocx clear
clear

/*
local folder "$resultsDir"

local files : dir "`folder'" files "*10*ss.emf"

// Create a new Word document
putdocx begin

// Loop through each .emf file and insert it into the Word document
foreach file of local files {
    
    // Insert the .emf graph into the Word document
    putdocx paragraph
    putdocx image "`folder'/`file'"
    
    // Optionally, you can add titles or additional text for each graph
    putdocx text ("Graph: `file'")
}

// Save and close the Word document
putdocx save "$resultsDir/stata_graphs_10.docx", replace

putdocx clear

local folder "$resultsDir"

local files : dir "`folder'" files "*small*ss.emf"

// Create a new Word document
putdocx begin

// Loop through each .emf file and insert it into the Word document
foreach file of local files {
    
    // Insert the .emf graph into the Word document
    putdocx paragraph
    putdocx image "`folder'/`file'"
    
    // Optionally, you can add titles or additional text for each graph
    putdocx text ("Graph: `file'")
}

// Save and close the Word document
putdocx save "$resultsDir/stata_graphs_smallconf.docx", replace

putdocx clear

local folder "$resultsDir"

local files : dir "`folder'" files "*no*ss.emf"

// Create a new Word document
putdocx begin

// Loop through each .emf file and insert it into the Word document
foreach file of local files {
    
    // Insert the .emf graph into the Word document
    putdocx paragraph
    putdocx image "`folder'/`file'"
    
    // Optionally, you can add titles or additional text for each graph
    putdocx text ("Graph: `file'")
}

// Save and close the Word document
putdocx save "$resultsDir/stata_graphs_noconf.docx", replace

