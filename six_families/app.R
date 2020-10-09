##production
library(stringr)
library(shiny)
library(ggplot2)
library(RColorBrewer)

#load distance matrix
dat = read.delim(paste0("data/mucin_TcYC6",".fasta.mat"), header = FALSE, sep = " ", quote = "\"", dec = ".", fill = TRUE)
dat2 <- dat[,-1]
rownames(dat2) <- dat[,1]
colnames(dat2) <- dat[,1]

#load loc
loc = read.delim(paste0("data/mucin_TcYC6",".fasta.mat.loc.tab"), col.names=c("gene_ID","chr","start","end","strand","annotation"), header = FALSE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
rownames(loc) <- loc[,1]
loc$length=loc$end-loc$start
loc$geneID_contig_length_annotation=paste(loc$gene_ID,loc$chr,loc$length,loc$annotation,sep=' ') 

#MDS
d <- as.dist(dat2)
mds.coor <- cmdscale(d)
colnames(mds.coor)=c("x","y")
plot_mat = merge(mds.coor,loc, by = "row.names" )
x_lims=c(min(plot_mat$x),max(plot_mat$x))
y_lims=c(min(plot_mat$y),max(plot_mat$y))

#shiny visualization

ui <- fluidPage(
  
  tags$head(tags$style('
     #my_tooltip {
      position: absolute;
      width: 600px;
      z-index: 100;
      padding: 0;
     }
  ')),
  
  tags$script('
    $(document).ready(function() {
      // id of the plot
      $("#distPlot").mousemove(function(e) { 

        // ID of uiOutput
        $("#my_tooltip").show();         
        $("#my_tooltip").css({             
          top: (e.pageY + 5) + "px",             
          left: (e.pageX + 5) + "px"         
        });     
      });     
    });
  '),
  fluidRow( column(width=4, selectInput("gFamily", "Select gene family:", selected="mucin_TcYC6", choices = c("DGF-1_TcBrA4","GP63_TcBrA4","MASP_TcBrA4","mucin_TcBrA4","RHS_TcBrA4","ts_TcBrA4","DGF-1_TcYC6","GP63_TcYC6","MASP_TcYC6","mucin_TcYC6","RHS_TcYC6","ts_TcYC6",
                                                                                                              "DGF-1_TcBrA4_fullAlnDist","GP63_TcBrA4_fullAlnDist","MASP_TcBrA4_fullAlnDist","mucin_TcBrA4_fullAlnDist","RHS_TcBrA4_fullAlnDist","ts_TcBrA4_fullAlnDist","DGF-1_TcYC6_fullAlnDist","GP63_TcYC6_fullAlnDist","MASP_TcYC6_fullAlnDist","mucin_TcYC6_fullAlnDist","RHS_TcYC6_fullAlnDist","ts_TcYC6_fullAlnDist",
                                                                                                              "DGF-1_TcBrA4_chrOnly_fullAlnDist","GP63_TcBrA4_chrOnly_fullAlnDist","MASP_TcBrA4_chrOnly_fullAlnDist","mucin_TcBrA4_chrOnly_fullAlnDist","RHS_TcBrA4_chrOnly_fullAlnDist","ts_TcBrA4_chrOnly_fullAlnDist","DGF-1_TcYC6_chrOnly_fullAlnDist","GP63_TcYC6_chrOnly_fullAlnDist","MASP_TcYC6_chrOnly_fullAlnDist","mucin_TcYC6_chrOnly_fullAlnDist","RHS_TcYC6_chrOnly_fullAlnDist","ts_TcYC6_chrOnly_fullAlnDist"
  ))),
            column(width=4, selectInput("var_y", "Mouse hover-over shows:",selected="geneID_contig_length_annotation", choices = c(names(plot_mat)[4:11],names(plot_mat)[2:3]))),
            #column(width=3, textInput('Contig', "Contig: (case sensitive)", value = "all", width = NULL, placeholder = "all")),
            #column(width=3, selectInput('Contig', "Contig:", selected = "all", choice=c("all",paste0(unique(loc$chr)))))
            uiOutput("chrSelection"),
            
  
  ),
  fluidRow(
            column(width=3, textInput('glen', 'gene length filter', value = "none", placeholder = "none")),
            column(width=3, textInput('dimX', 'Plot Width', value = "900", placeholder = 900)),
            column(width=3, textInput('dimY', 'Plot Height', value = "600", placeholder = 600))
            ),
  fluidRow(          
            column(width=10, textAreaInput("geneList", "Showing only genes in the list below: (one gene per line)", value = "all", width = 300,
                                          height = 200, cols = NULL, rows = NULL, placeholder = "all",
                                          resize = NULL)
                   ),
  ),            

  fluidRow(column(width=3, actionButton("applybutton", "Click to apply settings"))),
  
  fluidRow(
  plotOutput("distPlot", hover = hoverOpts(id="plot_hover",delay=0)),
  uiOutput("my_tooltip")
  )
)

server <- function(input, output) {

  plotHeight <- reactive(input$dimY)    
  plotWidth <- reactive(input$dimX) 
  
  output$distPlot <- renderPlot({
        isolate({
                ggplot(plot_mat, aes_string('x', 'y', color='chr',shape='chr')) + 
                  geom_point()+
                  scale_shape_manual(values=c(0:25,0:9,0:25,0:9,0:25,0:9,0:25,0:9,0:25,0:9))+
                  labs(title=paste0("MDS representation of pairwise genetic distances\n",input$gFamily),
                       x ="X coordinate in MDS space", y = "Y coordinate in MDS space")+
                  theme(plot.title = element_text(face="bold", size=15, hjust=0.5))+
                  theme(axis.title.x = element_text( size=12))+
                  theme(axis.title.x = element_text( size=12))
             
        })
       }, height= 600, width= 900
          
    )
  

  output$my_tooltip <- renderUI({
    hover <- input$plot_hover 
    y <- nearPoints(plot_mat, input$plot_hover, maxpoints=1)[input$var_y]
    req(nrow(y) != 0)
    verbatimTextOutput("vals")
  })
  
  output$vals <- renderPrint({
    hover <- input$plot_hover 
    y <- nearPoints(plot_mat, input$plot_hover)[input$var_y]
    req(nrow(y) != 0)
    y
  })    
  
  #populate chr selection
  output$chrSelection <- renderUI({
    column(width=3, selectInput('Contig', "Contig:", selected = "all", choice=c("all",paste0(unique(loc$chr)))))
  })
  
  observeEvent(input$applybutton,{
    isolate({  
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Progress: ", value = 0)
    progress$inc(1/5, detail = paste("loading distance matrix"))

      
      print(input$gFamily)
      print(input$Contig)
      print(input$geneList)
      gFamily_selection=input$gFamily
      Contig_selection=input$Contig
      geneList=input$geneList
      
      #load distance matrix
      dat = read.delim(paste0("data/",input$gFamily,".fasta.mat"), header = FALSE, sep = " ", quote = "\"", dec = ".", fill = TRUE)
      dat2 <- dat[,-1]
      rownames(dat2) <- dat[,1]
      colnames(dat2) <- dat[,1]
      
      progress$inc(1/5, detail = paste("loading gene annotation and applying various filters/settings"))

      
      #load loc
      loc = read.delim(paste0("data/",input$gFamily,".fasta.mat.loc.tab"), col.names=c("gene_ID","chr","start","end","strand","annotation"), header = FALSE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
      rownames(loc) <- loc[,1]
      loc$length=loc$end-loc$start
      loc$geneID_contig_length_annotation=paste(loc$gene_ID,loc$chr,loc$length,loc$annotation,sep=' ') 
 

      #length filter
      gene_lengt_CO = input$glen
      print(gene_lengt_CO)
      print(grepl("^[[:digit:]]",gene_lengt_CO))
      print(nrow(loc))
      if (grepl("^[[:digit:]]",gene_lengt_CO))
      {
        idx <- loc$length >= as.numeric(gene_lengt_CO)
        loc = loc[idx,]
        print(nrow(loc))
      }
      print(nrow(loc))
      
      #populate chr selection
      chr_selection = unique(loc$chr)
      output$chrSelection <- renderUI({
        column(width=3, selectInput('Contig', "Contig:", selected = "all", choice=c("all",paste0(chr_selection))))
      })
      
      #print(unique(loc$chr))
      
      #get species
      specie = strsplit(input$gFamily,"_")[[1]][2]
 
      #get axis limits
      d <- as.dist(dat2)
      mds.coor <- cmdscale(d)
      colnames(mds.coor)=c("x","y")
      plot_mat = merge(mds.coor,loc, by = "row.names" )
      x_lims=c(min(plot_mat$x),max(plot_mat$x))
      y_lims=c(min(plot_mat$y),max(plot_mat$y))
      
      progress$inc(1.5/5, detail = paste("rendering plot"))
      
      #MDS
      d <- as.dist(dat2)
      mds.coor <- cmdscale(d)
      colnames(mds.coor)=c("x","y")
      plot_mat = merge(mds.coor,loc, by = "row.names" )

      
      
      #filter contigs
      if (input$Contig != "all" && input$Contig != "All")
      {
        loc=loc[which(loc$chr==paste0(Contig_selection)),]
        plot_mat=plot_mat[which(plot_mat$chr==paste0(Contig_selection)),]
      }
      
      #filter genes

      if (geneList != "all" && geneList != "All" && geneList !="")
      {
        #parse input gene list
        geneList=gsub("[\t\n]+","\t", geneList)
        print(geneList)
        geneList=strsplit(geneList,"\t")
        print(geneList[[1]])
        
        #filter
        idx=match(geneList[[1]], plot_mat$gene_ID)
        loc=loc[idx,]
        plot_mat=plot_mat[idx,]        
      }
      

      progress$inc(1/5, detail = paste("rendering plot"))
      
      #print(plot_mat)   
      
      output$distPlot <- renderPlot({
        ggplot(plot_mat, aes_string('x', 'y', color='chr',shape='chr')) + 
          geom_point()+
          scale_shape_manual(values=c(0:25,0:9,0:25,0:9,0:25,0:9,0:25,0:9,0:25,0:9))+
          xlim(x_lims)+
          ylim(y_lims)+
          labs(title=paste0("MDS representation of pairwise genetic distances\n",gFamily_selection),
               x ="X coordinate in MDS space", y = "Y coordinate in MDS space")+
          theme(plot.title = element_text(face="bold", size=15, hjust=0.5))+
          theme(axis.title.x = element_text( size=12))+
          theme(axis.title.x = element_text( size=12))
      }, height= as.numeric(input$dimY), width= as.numeric(input$dimX))
      
      output$my_tooltip <- renderUI({
        hover <- input$plot_hover 
        y <- nearPoints(plot_mat, input$plot_hover, maxpoints=1)[input$var_y]
        req(nrow(y) != 0)
        verbatimTextOutput("vals")
      })
      
      output$vals <- renderPrint({
        hover <- input$plot_hover 
        y <- nearPoints(plot_mat, input$plot_hover)[input$var_y]
        req(nrow(y) != 0)
        y
      })  
    })
  })
}
shinyApp(ui = ui, server = server)

#library(rsconnect)
#deployApp()