## A Shiny app that visualizes clusters of gene family members based on genetic distance
### Multi-dimensional scaling is used to map genes onto a two-dimensional plot

###
This tool and associated data is reported in the following article:   
https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1009254

###
This tool is hosted at   
https://shiny.ctegd.uga.edu/  

### Instructions
Open one of the app.R files in Rstudio and hit the "run app" button

### Folder Structure:
```

.
├── six_families          # Mucins, MASP, ts, DGF-1, RHS, GP63 (k-tuple distance) 
├── ts                    # ts only, Brazil A4 and Y C6 (k-tuple distance)
├── ts_fullAlnDist        # ts only, Brazil A4 and Y C6 (full-alignment distance)
├── ts_fullAlnDist_v2     # ts only, Brazil A4 and Y C6 (refined full-alignment distance)
└── README.md

