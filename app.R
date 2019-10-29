

# Load libraries
library(tablerDash)
library(rols)
library(shinyEffects)
library(echarts4r)
library(shinyWidgets)
library(shinyjs)
library(karyoploteR)
library(DT)
library(XML)
library(gghighlight)
library(ReactomePA)
library(shinycssloaders)
library(plotly)
library(waiter)
library(DescTools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(DOSE)
library(enrichplot)
library(rentrez)
library(reactable)
library(ggridges)
library(UpSetR)
library(randomForest) # delete in case of using an alternative model
library(chromPlot)
import::from(Gviz, "IdeogramTrack")
import::from(Gviz, 'plotTracks')
library(tidyverse)
library(shiny)


load('local_data.RData')

source('functions.R')
# load('local_data.RData')
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
# file.remove('local_data.RData')
# save(hgcn_genes, df_enhancers, lncrna_coord, lncrna, tad, gtex, hpa, hpo_genes, vector_hp, vector_term, cnv_df,
#      select, cnv_df, model1, panel_total, denovo, file = "local_data.RData")



# Files needed to run the app

# hgcn_genes - List all the genes with scores
# df_enhancers - List enhancers with phast100vert  + phast46pla + phast46pri
# lncrna_coord - coordinates lncRNA
# lncrna - Features per lncRNA
# tad - List of TADs with coordinates
# gtex - gene expression data across tissues
# hpa - human protein expression
# hpo - Human phenotype ontology


theme_fancy <- function() {
  theme_minimal(base_family = "Asap Condensed") +
    theme(panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(size=22))
  
}



human_chrom <- hg19 <- list('chr1' = 1, 'chr2' = 2,'chr3' = 3,'chr4' = 4,'chr5' = 5,'chr6' = 6,'chr7' = 7,'chr8' = 8,'chr9' = 9,'chr10' = 10,'chr11' = 11,'chr12' = 12,'chr13' = 13,
                            'chr14' = 14,'chr15' = 15,'chr16' = 16,'chr17' = 17,'chr18' = 18,'chr19' = 19, 'chr20' = 20, 'chr21' = 21,  'chr22' = 22,
                            'chrX' = 'X','chrY' = 'Y')

hg_cytoBandIdeo <- chromPlot::hg_cytoBandIdeo




# app
shiny::shinyApp(
  
  
  
  ui = tablerDashPage(
    # enable_preloader = TRUE,
    # loading_duration = 4,
    

    
    
    navbar = tablerDashNav(
      
      navMenu = tablerNavMenu(

        tablerNavMenuItem(
          tabName = "overview",
          icon = "box",
          "Overview"
        ),
        tablerNavMenuItem(
          tabName = "fa",
          icon = "box",
          "Functional analysis"
        ),
        # tablerNavMenuItem(
        #   tabName = "model",
        #   icon = "box",
        #   "Model organism"
        # ),
        tablerNavMenuItem(
          tabName = "reg_region",
          icon = "box",
          "Regulatory regions"
        ),
        tablerNavMenuItem(
          tabName = "tissue",
          icon = "box",
          "Tissue-specificity"
        ),
        tablerNavMenuItem(
          tabName = "disease",
          icon = "box",
          "Disease-specificity"
        ),
        tablerNavMenuItem(
          tabName = "pubmed",
          icon = "book",
          "Biomedical literature"
        ),
        
        tablerNavMenuItem(
          tabName = "cnv_ngs",
          icon = "book",
          "CNVs NGS"
        ),
        tablerNavMenuItem(
          tabName = "down_report",
          icon = "download",
          "Automated report"
        )
        # tablerNavMenuItem(
        #   tabName = "Docu",
        #   icon = "book",
        #   "Documentation"
        # )
        
        
        ),
      id = "mymenu",
      src = "https://www.onlinelogomaker.com/applet_userdata/version2/5/0/18611424/projects/18611424.png",
      uiOutput('ref_user_filter_genes'),
      uiOutput('n_genes_added'),
      uiOutput('ref_user_genes'),
      uiOutput('ref_user_length')
      
      # tablerDropdown(
      #   tablerDropdownItem(
      #     title = "Item 1 title",
      #     href = "http://google.com",
      #     url = "https://image.flaticon.com/icons/svg/1301/1301804.svg",
      #     status = "danger",
      #     date = "now",
      #     "This is the first dropdown item"
      #   ),
      #   tablerDropdownItem(
      #     url = "https://image.flaticon.com/icons/svg/1301/1301809.svg",
      #     status = "warning",
      #     "This is the second dropdown item",
      #     date = "yesterday"
      #   ),
      #   tablerDropdownItem(
      #     title = "Item 3 title",
      #     "This is the third dropdown item"
      #   )
      # )
    ),
    footer = tablerDashFooter(
      tablerIcon(name = "maestro", lib = "payment"),
      tablerIcon(name = "mastercard", lib = "payment"),
      copyrights = "@David Granjon, 2019"
    ),
    title = "CNVxplorer - A web tool for the clinical interpretation of CNVs",
    body = tablerDashBody(
      
      tags$head(tags$script('
  $(document).on("shiny:sessioninitialized", function(event) {
    $(\'a[data-value="Page1"]\').tab("show");
  });
')),
      
      tablerTabItems(
        

        
        tablerTabItem(
          tabName = "overview",
          
          use_waiter(),
          setZoom(class = "card"),
          chooseSliderSkin("Nice"),
          
          fluidRow(
            column(
              width = 3,
              selectizeInput(inputId = 'input_chrom', label = 'Genomic interval - Chromosome', choices = human_chrom,
                             selected = NULL, multiple = FALSE,
                             options = NULL),
              prettyRadioButtons(
                inputId = "input_geno_karyo",
                label = "Choose:", 
                choices = c("Genomic coordinates", "G banding"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              uiOutput('choose_geno_karyo1'),
              uiOutput('choose_geno_karyo2'),
              # tags$hr(),
              actionBttn(
                inputId = "start_analysis",
                label = "Start!",
                color = "success",
                style = "material-flat",
                # icon = icon("sliders"),
                block = TRUE
              ),
              tags$br(),
              tags$br(),
              tags$br()
            ),
            column(6,
                   plotOutput('plot_chrom', height = 200)),
            column(
              width = 3,

              # prettyRadioButtons(
              #   inputId = "snv_yes_no",
              #   label = "Upload SNVs file (.bed)?",
              #   choices = c("No", "Yes"),
              #   inline = TRUE,
              #   status = "primary",
              #   fill = TRUE
              # ),
              # conditionalPanel(
              #   condition = "input.snv_yes_no == 'Yes'",
              #   fileInput("file_snv", label = h5("Upload file:"))
              # ),
              # hr(),
              # uiOutput('n_snv'),
              
              uiOutput('n_genes'),
              uiOutput('score_rf'),
              
              uiOutput("info")
            )
          ),
          fluidRow(
            # column(
            #   width = 6,
            #   tablerCard(
            #     title = "Plots",
            #     zoomable = FALSE,
            #     closable = FALSE,
            #     options = tagList(
            #       switchInput(
            #         inputId = "enable_distPlot",
            #         label = "Plot?",
            #         value = TRUE,
            #         onStatus = "success",
            #         offStatus = "danger"
            #       )
            #     ),
            #     plotOutput("distPlot"),
            #     status = "info",
            #     statusSide = "left",
            #     width = 12,
            #     footer = tagList(
            #       column(
            #         width = 12,
            #         align = "center",
            #         sliderInput(
            #           "obs",
            #           "Number of observations:",
            #           min = 0,
            #           max = 1000,
            #           value = 500
            #         )
            #       )
            #     )
            #   )
            # ),
            
            tablerCard(
              title = "Gene panel",
              DTOutput("dgenes"),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE,
              options = tagList(
                # actionBttn(
                #   inputId = "Id111",
                #   label = "Reset table", 
                #   style = "minimal",
                #   color = "primary"
                # ),
                downloadButton("download_dgenes", "Download table")
                
                
            )
            ),
            tablerCard(
              title = "Comparison scores",
              collapsible = FALSE,
              closable = FALSE,
              
              plotlyOutput("dot_comparison"),
              width = 12,
              overflow = TRUE
            ),
            column(6,
            tablerCard(
              title = "Clinical panels (Genomics England PanelApp)",
              DTOutput("df_fisher"),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE,
              options = tagList(
                # actionBttn(
                #   inputId = "Id111",
                #   label = "Reset table", 
                #   style = "minimal",
                #   color = "primary"
                # ),
              )
            ),
            

            
            tablerCard(
              title = "Gene(s) selected",
              # DTOutput("df_fisher_selected"),
              textOutput("df_fisher_selected"),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE

              
            )),
            
            tablerCard(
              title = "Gene-level databases",
              DTOutput("df_fisher_two"),
              width = 6,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE,
              options = tagList(
              )
            ),
            
            column(width = 6,
                   
                   uiOutput('n_cnv_patho'),
                   uiOutput('n_cnv_nopatho'),
                   uiOutput('n_clinvar'),
                   uiOutput('n_gwas')),
            column(width = 6,
                   
                   uiOutput('n_dev'),
                   uiOutput('n_omim'),
                   # uiOutput('n_pubmed'),
                   uiOutput('n_pli'),
                   uiOutput('n_hi'))
            
            
          ),
          fluidRow(
            # tablerCard(
            #   title = "Heatmap - Genome-wide percentile",
            #   plotlyOutput('heatmap'),
            #   width = 12,
            #   overflow = TRUE
            # ),
            # tablerCard(
            #   title = "Overlapping of genes and CNV",
            #   plotOutput('plotp_overlap'),
            #   width = 6,
            #   overflow = TRUE
            # ),
            tablerCard(
              title = "De novo variants (denovo-db)",
              DTOutput('df_de_novo'),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE,
              options = tagList(
                downloadButton("download_de_novos", "Download table")
                
              )
            ),
            tablerCard(
              title = "Comparison CNV size with other CNVs databases (gnomAD, DGV and DECIPHER)",
              plotOutput('plot_size'),
              width = 12,
              overflow = TRUE,
              collapsible = FALSE,
              closable = FALSE,       
              options = tagList(

                pickerInput(
                  inputId = "select_density",
                  width = 'fit',
                  # label = "sdasda", 
                  choices = c("global", "local")
                )
                
              )
            ),
            tablerCard(
              title = "Overlapping with pathogenic CNVs (DECIPHER)",
              DTOutput('df_overlap_cnvs'),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE,
              options = tagList(
                downloadButton("download_df_overlap_cnvs", "Download table")
                
              )
            ),
            tablerCard(
              title = "Overlapping with non-pathogenic CNVs (gnomAD, DGV)",
              DTOutput('df_overlap_cnvs_nonpatho'),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE,
              options = tagList(
                downloadButton("download_df_overlap_cnvs_nonpatho", "Download table")
                
              )
            ),
            tablerCard(
              title = "Barplot ",
              plotOutput('plot_cnv_bar'),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE
            ),
            # column(9,
            tablerCard(
              title = "Intersection - CNV and CNV pathogenics",
              DTOutput('df_intersection'),
              width = 5,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE,
              options = tagList(
                uiOutput('ui_intersect')
              )
            ),
            # useShinyjs()
            # ),
            tablerCard(
              title = "Plot intersection",
              plotOutput('plot_intersection'),
              width = 12,
              collapsible = FALSE,
              closable = FALSE,            
              overflow = TRUE
            )
            # tablerCard(
            #   title = "Comparison CNV size with other CNVs (gnomAD, DGV and DECIPHER)",
            #   DTOutput('test_reading_snv_file'),
            #   width = 12,
            #   overflow = TRUE
            # )

          )
          # tablerCard(
          #   title = "Functional analysis",
          #   width = 12,
          #   zoomable = FALSE,
          #   closable = FALSE,
          #   plotOutput('func_an1alysis') %>% withSpinner(type = 5),
          #   options = tagList(
          #     switchInput(
          #       inputId = "enable1_func_analysis",
          #       label = "Run?",
          #       value = FALSE,
          #       onStatus = "success",
          #       offStatus = "danger"
          # 
          #     ),
          #     pickerInput(
          #       inputId = "sign_1vline",
          #       label = tags$b("P.value threshold:"),
          #       choices = c("0.05", "0.01", "0.005"),
          #       width = 130),
          #     prettyRadioButtons(
          #       inputId = "choos1e_go",
          #       label = tags$b("Select one option:"),
          #       choices = c("biological process", "molecular function", "cellular component"),
          #       inline = TRUE,
          #       status = "primary",
          #       fill = TRUE
          #     )
          #   )),
          # tablerCard(
          #   title = "Gene-disease association",
          #   plotOutput('func_anal1ysis_diseases') %>% withSpinner(type = 5),
          #   width = 12,
          #   zoomable = FALSE,
          #   closable = FALSE
          # 
          # 
          # )
          # tablerCard(
          #   title = "Gene-disease association",
          #   plotOutput('func_2nalysis_mesh') %>% withSpinner(type = 5),
          #   width = 12,
          #   zoomable = FALSE,
          #   closable = FALSE
          #   
          #   
          # )
        ),
        
        tablerTabItem(
          tabName = "fa",
          # fluidRow(
          #   
          #   tablerCard(
          #     
          #     title = 'P.value threshold',
          #     width = 3,
          #     pickerInput(
          #       inputId = "sign_vline",
          #       label = tags$b("P.value threshold:"), 
          #       choices = c("0.05", "0.01", "0.005"),
          #       width = 12)
          #     
          #     
          #   )
          #   
          #   
          # ),
          tablerCard(
            title = "Functional Profile",
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            uiOutput('plot_df') ,
            options = tagList(
            # column(width = 1
            # 
            #   # switchInput(
            #   #   inputId = "user_df_plot",
            #   #   label = "Table?",
            #   #   value = FALSE,
            #   #   onStatus = "success",
            #   #   offStatus = "danger"
            #   # )
            #   ),

              pickerInput(
                inputId = "user_level",
                label = tags$b("Level:"), 
                choices = 1:52,
                selected = 2,
                width = 130,
                options = list(
                  size = 10,
                  `live-search` = TRUE)),
            tags$hr(),
              prettyRadioButtons(
                inputId = "choose_group_go",
                label = tags$b("Select one option:"), 
                choices = c("biological process", "molecular function", "cellular component"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
            prettyRadioButtons(
              inputId = "user_df_plot",
              label = tags$b("Select one option:"), 
              choices = c("Bar plot", "Table"),
              inline = TRUE, 
              status = "primary",
              fill = TRUE
            ),
            
            switchInput(
              inputId = "enable_group_go",
              label = "Run?",
              value = FALSE,
              onStatus = "success",
              offStatus = "danger",
              width = 'auto',
              size = 'mini'
              
            )
              
              
            )
          ),
          # tablerCard(
          #   title = "Gene ontology",
          #   width = 12,
          #   zoomable = FALSE,
          #   closable = FALSE,
          #   DTOutput('df_22go') %>% withSpinner(type = 5)
          # ),
          tablerCard(
            title = "Gene Ontology",
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            uiOutput('plot_enrichgo') %>% withSpinner(type = 5),
            options = tagList(
              pickerInput(
                inputId = "sign_vline",
                label = tags$b("P.value threshold:"), 
                choices = c("0.05", "0.01", "0.005"),
                width = 130),
              prettyRadioButtons(
                inputId = "choose_go",
                label = tags$b("Select one option:"), 
                choices = c("biological process", "molecular function", "cellular component"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              prettyRadioButtons(
                inputId = "table_plot_go",
                label = tags$b("Select one option:"), 
                choices = c("Bar plot", "Table"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              switchInput(
                size = 'mini',
                inputId = "enable_func_analysis",
                label = "Run?",
                value = FALSE,
                onStatus = "success",
                offStatus = "danger",
                width = 'auto'
                
              )
            )),
          tablerCard(
            title = "Pathway analysis",
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            uiOutput("ui_path") %>% withSpinner(type = 5),
            options = tagList(
              pickerInput(
                inputId = "sign_vline_path",
                label = tags$b("P.value threshold:"), 
                choices = c("0.05", "0.01", "0.005"),
                width = 130),
              prettyRadioButtons(
                inputId = "kegg_reactome",
                label = tags$b("Select a database:"), 
                choices = c("Reactome", "KEGG"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              prettyRadioButtons(
                inputId = "table_path_go",
                label = tags$b("Select one option:"), 
                choices = c("Bar plot", "Table"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              switchInput(
                size = 'mini',
                inputId = "enable_path_analysis",
                label = "Run?",
                value = FALSE,
                onStatus = "success",
                offStatus = "danger",
                width = 'auto'
                
              )
          )),
          tablerCard(
            title = "Gene-disease association",
            uiOutput('ui_do') %>% withSpinner(type = 5),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            options = tagList(
              pickerInput(
                inputId = "pvalue_do",
                label = tags$b("P.value threshold:"), 
                choices = c("0.05", "0.01", "0.005"),
                width = 130),
              prettyRadioButtons(
                inputId = "table_plot_do",
                label = tags$b("Select one option:"), 
                choices = c("Gene-disease plot", "Table"),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              switchInput(
                size = 'mini',
                inputId = "enable_do_analysis",
                label = "Run?",
                value = FALSE,
                onStatus = "success",
                offStatus = "danger",
                width = 'auto'
                
              )
            )
          )
          
        ),
        tablerTabItem(
          tabName = "reg_region",
          fluidRow(
            # tablerCard(title = 'Select a region:',
            #            uiOutput('gen2e_2tissue'),
            #            width = 3),
            uiOutput('n_enhancer_total'),
            uiOutput('n_tads'),
            uiOutput('n_lncrna'),
            
            # uiOutput('n_enhancer_inside'),
            uiOutput('n_enhancer'),
            
            # uiOutput('redund_n_enhancer'),
            # tablerCard(title = 'Include target-genes enhancers',
            #            uiOutput('switch_enhancers'),
            #            collapsible = FALSE,
            #            closable = FALSE,
            #            width = 3),
            # uiOutput('switch_enhancers'),
            
            
            tablerCard(title = 'Enhancers disrupted',
                       collapsible = FALSE,
                       closable = FALSE,
                       dataTableOutput('df_enhancer'),
                       width = 12,
                       options = tagList(
                        uiOutput('n_filtered_enhancers'),
                        uiOutput('switch_enhancers')
                        
                         
                         
                         ))),
          tablerCard(title = 'Phast100way histogram',
                     plotOutput('histo_enhancer'),
                     collapsible = FALSE,
                     closable = FALSE,
                     width = 6),
          tablerCard(title = 'Long noncoding RNA (lncRNA) disrupted',
                     DTOutput('lncrna_df'),
                     collapsible = FALSE,
                     closable = FALSE,
                     width = 12),
          tablerCard(title = 'Topologically Associating Domains (TADs) disrupted',
                     DTOutput('df_tads'),
                     collapsible = FALSE,
                     closable = FALSE,
                     width = 12)
          
        ),
        
        tablerTabItem(
          tabName = "cnv_ngs",
          fluidRow(
            # column(3,
            # uiOutput('n_upload_cnv'),
            # uiOutput('n_upload_filter_cnv')),
            column(4,
            tablerCard(title = 'Input file',
                       width = 12,
                       collapsible = FALSE,
                       closable = FALSE,
            fileInput("upload_bed_file", label = h5(strong("CNVs regions (.bed file)")), 
                      accept = c(".bed"))),
            tablerCard(title = 'Filter by size',
                       width = 12,
                       collapsible = FALSE,
                       closable = FALSE,
                       sliderTextInput(
                         inputId = "filter_length",
                         label = "",
                         choices = c(1, 10, 100, 500, 1000, 10000),
                         grid = TRUE
                       )),
            tablerCard(title = 'Filter by pathogenicity score',
                       width = 12,
                       collapsible = FALSE,
                       closable = FALSE,
                       sliderTextInput(
                         inputId = "filters_length",
                         label = "",
                         choices = c(1, 10, 100, 500, 1000, 10000),
                         grid = TRUE
                       )),
            
            tablerCard(title = 'Exclude benign CNVs',
                       width = 12,
                       collapsible = FALSE,
                       closable = FALSE,
                       sliderTextInput(
                         inputId = "filters_laength",
                         label = "",
                         choices = c(1, 10, 100, 500, 1000, 10000),
                         grid = TRUE
                       )),
            tablerCard(title = 'Filter by number of genes ',
                       width = 12,
                       collapsible = FALSE,
                       closable = FALSE,
                       sliderTextInput(
                         inputId = "filters_leangth",
                         label = "",
                         choices = c(1, 10, 100, 500, 1000, 10000),
                         grid = TRUE
                       ))

            ),

            tablerCard(title = 'List of CNVs',
                       collapsible = FALSE,
                       closable = FALSE,
                       dataTableOutput('df_upload'),
                       width = 8
                       ))

          
        ),
        
        tablerTabItem(
          tabName = "gene",
          fluidRow(
            tablerCard(title = '19,281 protein-coding genes ',
                       reactableOutput("table_all_genes"),
                       width = 12)),
          tablerCard(title = 'RNA Expression (GTEx)',
                     DTOutput('lncrn2a_df'),
                     width = 12)
          
        ),
        
        tablerTabItem(
          tabName = "tissue",
          fluidRow(
            tablerCard(title = 'Select a gene:',
                       uiOutput('gene_tissue'),
                       width = 3),
            tablerCard(title = 'Protein Expression (Human Protein Atlas)',
                       DTOutput('tissue_hpa'),
                       width = 9)),
          tablerCard(title = 'RNA Expression (GTEx)',
                     plotlyOutput('tissue_gtex'),
                     width = 12)
          
        ),
        tablerTabItem(
          tabName = "disease",
          fluidRow(
            tablerCard(title = 'Phenotype terms',
            width = 9,
            multiInput(
              inputId = "chosen_hp",
              label = "", 
              choices = NULL,
              width = '100%',
              choiceNames = vector_term,
              choiceValues = vector_hp
              
            )),
            column(width = 3,
            uiOutput('n_hp_chosen'),
            uiOutput('check_genes_hp')),
            tablerCard(title = 'Genes associated with the phenotype terms',
                       DTOutput('df_check_hp_genes'),
                       width = 12),
            tablerCard(title = 'Intersection of phenotype terms',
                       plotOutput('plot_upset'),
                       width = 12),
            # width = 6
            # ),
            # tablerCard(title = 'Select a gene:',
            #            uiOutput('n_pub2med'),
            #            width = 3),
            tablerCard(title = 'Pubmed articles associated with the region',
                       collapsible = FALSE,
                       closable = FALSE,
                       DTOutput('disease2_pubmed'),
                       width = 12,
                       options = tagList(
                         downloadButton("download_pubmed", "Download table")
                       )))
        ),
        tablerTabItem(
          tabName = "model",
          tablerCard(title = 'Phenotypes associated with the list of genes',
                     echarts4rOutput('agg_model'),
                     width = 12),
          fluidRow(
            # tablerCard(title = 'Select a gene:',
            #            uiOutput('gene_model'),
            #            width = 3),
            tablerCard(title = 'Genes associated with phenotypes in mouse (MGI)',
                       DTOutput('model_genes'),
                       width = 12))
          
        ),
        tablerTabItem(
          tabName = "pubmed",
          fluidRow(
            tablerCard(
              title = "Pubmed articles associated with the region",
              
              collapsible = FALSE,
              closable = FALSE,
              DTOutput('disease_pubmed'),
              overflow = TRUE,
              width = 12
            ),
            tablerCard(
              title = "Abstract",
              
              DTOutput("abstract_pubmed"),
              width = 12,
              overflow = TRUE
            )
            # tablerCard(title = 'Download report (.html)',
            #            downloadButton('button55_download', label = "Download"),                       
            #            width = 6),
            
            # tablerCard(title = '',
            #            DTOutput('mode2l55_genes'),
            #            width = 6)
            
            )
          
        ),
        tablerTabItem(
          tabName = "down_report",
          fluidRow(
            tablerCard(
              title = "Personalize your report:",
              
              # multiInput(
              #   inputId = "Id010",
              #   label = "Countries :", 
              #   choices = NULL,
              #   choiceNames = c('Name clinician', 'Age', 'Sex', 'Add comment'), 
              #   choiceValues = as.character(1:4)
              #   ),
              fluidRow(
              column(width = 2,
                     tags$b('Patient information'),
                     tags$hr(),
              # textInput("name_report", "Name:", "Doctor Requena"),
              # textInput("age_report", "Age:", "4"),
              dateInput("age_report", "Date of birth (mm/dd/yy):", value = "", format = "mm/dd/yy"),
              pickerInput("sex_report", "Sex:", c("Male", "Female"))),
              column(width = 2,
                     tags$b('Clinician information'),
                     tags$hr(),
                     textInput("name_report", "Name:", placeholder = '')),
                     # textInput("age_report", "Age:", "4"),
                     # dateInput("age_report", "Date of birth:", value = "2012-02-29", format = "mm/dd/yy"),
                     # pickerInput("sex_report", "Sex:", c("Male", "Female"))),
              column(width = 8,
              textAreaInput("comment_report", "Notes:", placeholder =  "Notes...", 
                            # width = '700px', 
                            heigh = '250px'))),
              overflow = TRUE,
              width = 7
            ),
            tablerCard(
              title = "Funnel overview",
              
              echarts4rOutput("funnel_genes"),
              width = 5,
              overflow = TRUE
            ),
            tablerCard(title = 'Download report (.html)',
                       downloadButton('button_download', label = "Download"),                       
                       width = 7)
           
          )
          
        )
        
        
        
        
        
        
        
        
        
      )
      
    )
  ),
  server = function(input, output, session) {
    
    
    
    
    gene_selected <- reactive({
      
      test4 <<- input$dgenes_rows_selected
      
      test_1993 <<- input$dgenes_state
      
      
      print(test4)
    })
    
    # observeEvent(input$start_analysis, {
    # 
    #   shinyjs::reset("start_analysis")
    # 
    # })
    
    # observe({
    #   input$start_analysis
    #   updateActionButton(session, 'start_analysis', value = 0)
    #   
    #   
    # })
    
    

    
    observeEvent(input$take_intersect, {
      
      
        df_tmp <- intersection_running() %>%
          select(id, start, end) %>%
          slice(input$df_intersection_rows_selected)
        
        test32134124122121421412 <<- df_tmp
        
        new_start <- df_tmp %>% pull(start)
        new_end <- df_tmp %>% pull(end)
      
      
      updateNumericInput(session, "int_start", value = new_start)
      updateNumericInput(session, "int_end", value = new_end)

    })
    
    
    coord_user <- eventReactive(input$start_analysis, {
      
      # req(input$start_analysis > 0)
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        coord_start <- input$int_start
        coord_end <-  input$int_end
        coord_chrom <- input$input_chrom
        
        coord_start <- as.numeric(coord_start)
        coord_end <- as.numeric(coord_end)
        
        
      } else {
        
        df_tmp <- chromPlot::hg_cytoBandIdeo %>%
          filter(Chrom == coord_user()[3]) %>%
          filter(Name == input$input_karyotype)
        
        coord_start <- df_tmp %>% select(Start) %>% pull()
        coord_end <-  df_tmp %>% select(End) %>% pull()
        coord_chrom <- input$input_chrom
        
        coord_start <- as.numeric(coord_start)
        coord_end <- as.numeric(coord_end)
    
      }
      
      # if (!is.null(input$take_intersect)) {
      #   if (input$take_intersect > 0) {
      #     
      
      # if (exists('input$df_intersection_rows_selected')) {
      #   if (base::exists("input$take_intersect")) {
      #     if (!is.null(input$df_intersection_rows_selected)) {
      # 
      #   test655 <<- input$take_intersect
      #   test92313122131321 <<- input$df_intersection_rows_selected
      #   
      #   df_tmp <- intersection_running() %>%
      #     select(id, start, end) %>%
      #     slice(input$df_intersection_rows_selected)
      #   
      # 
      #   coord_start <- df_tmp %>% pull(start)
      #   coord_end <-  df_tmp %>% pull(end)
      #   coord_chrom <- input$input_chrom
      # 
      #   }
      # }
        
      coord_start <- as.numeric(coord_start)
      coord_end <- as.numeric(coord_end)
      
      c_output <- c(coord_start, coord_end, coord_chrom)
      

      # shinyjs::reset('start_analysis')
      
    })
    
    
    
    
    output$gene_tissue <- renderUI({
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      input_data <- df_genes %>% select(gene) %>% pull()
      
      pickerInput(
        inputId = "input_gene_tissue",
        # label = "Select gene:", 
        choices = input_data,
        options = list(
          size = 5,
          `live-search` = TRUE))
      
      
      
    })
    
    
    output$gene_model <- renderUI({
      
      
      input_data <- data_selected() %>% select(gene) %>% pull()
      
      pickerInput(
        inputId = "input_gene_model",
        # label = "Select gene:", 
        choices = input_data,
        options = list(
          size = 10,
          `live-search` = TRUE))
      
      
      
    })
    
    
    
    check_cnv_df <- reactive({
      
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      df_output <- cnv_df %>%
        filter(chrom == chrom_coordinates) %>%
        rowwise() %>%
        # check this because i dont get why did i use two overlaps
        mutate(keep = c(start, end) %overlaps% c(start_coordinates, end_coordinates)) %>%
        filter(keep == TRUE) %>%
        select(-keep)
      
      test312 <<- df_output
      df_output
      
    })
    
    
    check_hp_genes <- reactive({
      
      req(input$start_analysis > 0)
      
      go <- Ontology("hp")
      
      
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      
      validate(
        need(!is.null(input$chosen_hp), "Please, select phenotype terms.")
      )
      
      test322 <<- input$chosen_hp
      
      hp_chosen <- input$chosen_hp
      genes_chosen <- df_genes %>% select(entrez_id) %>% pull()
      
      
      df_tmp <- hpo_genes %>%
        filter(hp %in% hp_chosen) %>%
        filter(entrez_id %in% genes_chosen)
      
      test9312321 <<- df_tmp
      
      df_tmp2 <- df_tmp %>%
        select(hp) %>%
        distinct() %>%
        mutate(description = map(hp, function(x) ifelse(length(term(go, x)@description) == 0, NA, term(go, x)@description))) 
        
      
      df_tmp <- df_tmp %>% 
        left_join(df_tmp2, by = 'hp')
      
      df_tmp
      
    })
    
    output$check_genes_hp <- renderUI({
      
      req(input$start_analysis > 0)
      
      test66 <<- check_hp_genes()
      
      tablerStatCard(
        value =  nrow(check_hp_genes()),
        title = "Number of genes associated with the phenotype(s)",
        # trend = -10,
        width = 12
      )
      
    })
    
    output$df_check_hp_genes <- renderDataTable({
      
      req(input$start_analysis > 0)
      
      test766 <<- check_hp_genes()
      
      validate(
        need(nrow(check_hp_genes()) != 0, "Not genes found.")
      )

      datatable(check_hp_genes(), options = list(
        pageLength = 5))
    })
    
    output$plot_upset <- renderPlot({
      
      test771 <<- check_hp_genes()
    get_upset(check_hp_genes())
      
      
    })
    
    output$n_hp_chosen <- renderUI({
      
      req(input$start_analysis > 0)
      
      tablerStatCard(
        value =  length(input$chosen_hp),
        title = "Clinical features selected",
        # trend = -10,
        width = 12
      )
      
    })
    
    output$n_cnv_patho <- renderUI({
      
      req(input$start_analysis > 0)
      
      data_tmp <- check_cnv_df() %>% filter(source == 'decipher') %>% nrow()
      
      tablerStatCard(
        value =  data_tmp,
        title = "Number of pathogenic CNVs",
        # trend = -10,
        width = 12
      )
      
      
    })
    
    output$n_cnv_nopatho <- renderUI({
      
      req(input$start_analysis > 0)
      
      data_tmp <- check_cnv_df() %>% filter(source != 'decipher') %>% nrow()
      
      tablerStatCard(
        value =  data_tmp,
        title = "Number of nonpathogenic CNVs",
        # trend = -10,
        width = 12
      )
      
      
    })
    
    
    query_pubmed <- reactive({
      
      req(input$start_analysis > 0)
      
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        
        start_coordinates <- coord_user()[1]
        end_coordinates <- coord_user()[2]
        chrom_coordinates <- coord_user()[3]
        
        query_region <-  chromPlot::hg_cytoBandIdeo %>%
          filter(Chrom %in% chrom_coordinates) %>%
          mutate(keep = map2_chr(Start, End, function(x,y) c(start_coordinates, end_coordinates) %overlaps% c(x,y))) %>%
          filter(keep == TRUE) %>%
          select(Name) %>%
          pull() %>%
          map_chr(function(x) paste0(chrom_coordinates, x)) %>%
          paste0(collapse = ' OR ')
        
      } else {
        query_region <- paste0(chrom_coordinates, input$input_karyotype)
        
      }
      
      test67 <<- query_region
      
      query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 200 )
      
      
      
    })
    
    output$n_pubmed <- renderUI({
      
      tablerStatCard(
        value =   length(query_pubmed()[['ids']]),
        title = "Number of articles found in Pubmed",
        width = 12
      )
      
    })
    
    
    running_pubmed <- reactive({
      
      
      
      # validate(
      #   need(length(query_pubmed()) != 0, "Please, select a gene in the datatable.")
      # )
      
      test21 <<- query_pubmed()
      query_tmp <- entrez_summary(db="pubmed", id= query_pubmed()[['ids']])
      
      test987 <<-query_tmp
      
      title <- unname(map_chr(query_tmp, function(x) x[["title"]]))
      n_cites <- unname(map_chr(query_tmp, function(x) x[["pmcrefcount"]]))
      first_author <- unname(map_chr(query_tmp, function(x) x[["sortfirstauthor"]]))
      last_author <- unname(map_chr(query_tmp, function(x) x[["lastauthor"]]))
      date_release <- unname(map_chr(query_tmp, function(x) x[["pubdate"]]))
      pmid_id <- unname(map_chr(query_tmp, function(x) x[["uid"]]))
      journal_id <- unname(map_chr(query_tmp, function(x) x[["fulljournalname"]]))
      
      
      df_output <- tibble(title = title, first_author = first_author, last_author = last_author, 
                          n_cites = n_cites, journal = journal_id, date_release = date_release, 
                          pmid = pmid_id)
      
      df_output <- df_output %>% mutate(n_cites = as.integer(ifelse(n_cites == '', 0, n_cites)))
      df_output
      
      
    })
    
    output$disease_pubmed <- renderDataTable({
      
      
      datatable(running_pubmed(), rownames = FALSE, filter = 'top', selection = 'single',

                colnames = c('Title','First author', 'Last author', 'N°cites','Journal', 'Published date', 'PMID' ),
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                  selection = 'single'
                  # columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
      
    })
    
    
    output$abstract_pubmed <- renderDataTable({
      
      validate(
        need(input$disease_pubmed_rows_selected != '', "Please, select an article from the datatable.")
      )
      
      test931 <<- input$disease_pubmed_rows_selected
      paper_selected <- running_pubmed() %>% slice(input$disease_pubmed_rows_selected)
      
      test45 <<- paper_selected
      
      fetch.pubmed <- entrez_fetch(db = "pubmed", id = test45 %>% pull(pmid), rettype = "xml", parsed = T)
      # Extract the Abstracts for the respective IDS.  
      abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
        xmlValue(xmlChildren(x)$Abstract))
      
      tmp_df <- tibble(title = paper_selected$title, abstract = abstracts[[1]]) %>% gather('Category', 'Info')
      tmp888 <<- tmp_df
      
      datatable(tmp_df, rownames = FALSE,
                options = list(dom = 't'))
      
    })
    
    
    
    data_selected_prev <- reactive({
      
      # req(input$start_analysis > 0)
      # req(input$start_analysis)
      req(coord_user())
      
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      test999 <<- coord_user()[1]
      test1000 <<- coord_user()[2]
      
      
      
      data_raw <- hgcn_genes %>% filter(chrom == chrom_coordinates)
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        data_raw <- data_raw  %>% mutate(keep = NA) %>%
          rowwise() %>%
          mutate(keep = c(start_position, end_position) %overlaps% c(start_coordinates, end_coordinates)) %>%
          filter(keep == TRUE) %>% 
          select(-keep) %>%
          ungroup()
        
        test0 <<- data_raw
        
      } else {
        
        data_raw <- data_raw  %>% mutate(keep = NA) %>%
          rowwise() %>%
          mutate(keep = c(start_position, end_position) %overlaps% c(start_coordinates, end_coordinates)) %>%
          filter(keep == TRUE) %>% 
          select(-keep) %>%
          ungroup()
      }
      # if (input$snv_yes_no == 'Yes') {
      #   
      #   snv_df <- reading_snv_file()
      #   
      # 
      #   for (i in 1:nrow(snv_df)) {
      #     
      #     df_add <- hgcn_genes %>% 
      #       filter(chrom == snv_df$chrom[i]) %>%
      #       mutate(keep = NA) %>%
      #       rowwise() %>%
      #       mutate(keep = c(start_position, end_position) %overlaps% c(snv_df$start[i], snv_df$end[i])) %>%
      #       filter(keep == TRUE) %>% 
      #       select(-keep) %>%
      #       ungroup()
      #     # check with snvs that overlap in more than 1 gene
      #     if (!df_add$gene %in% data_raw$gene) {
      #       data_raw <-  bind_rows(data_raw, df_add)
      #     }
      #   }}

      data_raw <- data_raw %>% 
        mutate(coordinates = paste0(chrom,':', start_position,'-', end_position)) %>% 
        mutate(oe = paste0(oe_lof, ' (', oe_lof_lower, '-',oe_lof_upper, ')')) %>%
        select(-transcript, -oe_lof, -oe_lof_lower, -oe_lof_upper, -vg, -ensembl_gene_id)
      
      test1521 <<- data_raw

      data_raw <- get_perc_overlap(data_raw, coord_user()[1], coord_user()[2])
      
      data_raw <- data_raw %>% distinct()
      
      test007 <<- data_raw
      
      data_raw

      return(data_raw)
      
    })
    
    data_selected_enhancers <- reactive({
      
      if (!is.null(input$enhancers_on_off)) {
        if (input$enhancers_on_off) {
      
        # Genes  mapped in CNV
        genes_cnv <- data_selected_prev()
        genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # Genes NOT mapped in CNV
        if (is.null(input$df_enhancer_rows_all)) {
          enhancers_df <- prev_enhancer()
        } else {
          enhancers_df <- prev_enhancer()[input$df_enhancer_rows_all,]
        }
        genes_no_cnv <- enhancers_df %>% select(gene) %>% distinct() %>% pull()
        genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        # ADAPT IT WHEN ADDING OMIM OR OTHERS!!!
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-oe_lof, -oe_lof_lower, -oe_lof_upper, -vg, -transcript, -ensembl_gene_id)

        } else {
          table_output <- tibble()
        }
      } else {
        table_output <- tibble()
      }
      table_output
    })
    
    
    data_selected <- reactive({
      
      
      table_output <- data_selected_prev() %>% 
        mutate(source = 'CNV') %>% 
        bind_rows(data_selected_enhancers() %>% mutate(source = 'Enhancer')) %>%
        mutate(source = as.factor(source))
      
      test2019 <<- table_output
      table_output
      
    })
    

    # data_selected <- reactive({
    #   
    #   
    #   if (is.null(input$dgenes_rows_all)) {
    #     df_output <- data_selected_prev()
    #   } else {
    #     df_output <- data_selected_prev()[input$dgenes_rows_all,]
    #   }
    # 
    # 
    #   # if (length(input$dgenes_rows_all) != nrow(data_selected_prev())) {
    #   # 
    #   #   data_selected <- data_selected_prev()
    #   #   data_selected <- data_selected[input$dgenes_rows_all,]
    #   # }
    # 
    # 
    #   df_output
    # 
    # 
    # })
    
    output$dgenes <- renderDataTable({
      
      
      
      server <- TRUE
      data_input <- data_selected() %>% select(-start_position, -end_position, -chrom)
      
      datatable(data_input, rownames = FALSE, filter = 'top', 
                selection = 'single',
                # extensions = 'Responsive',
                options = list(
                  pageLength = 5, 
                  # autoWidth = TRUE, 
                  style = 'bootstrap', 
                  list(searchHighlight = TRUE),
                  stateSave = FALSE
                  # colnames = c('Entrez id', 'Band', 'Gene', 'pLI', 'Database', 'CNV size', 'Percentage Overlap (%)')
                  # columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
        # formatStyle(
        #   'pLI',
        #   background = styleColorBar(c(0,1), '#ca7171'),
        #   backgroundSize = '100% 90%',
        #   backgroundRepeat = 'no-repeat',
        #   backgroundPosition = 'center'
        # )
    })
    
    output$score_references <- renderDataTable({
      
      
      test1945 <<- input$dgenes_rows_all
      
      
      datatable(ref_scores, rownames = FALSE,
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap',
                  selection = 'single',
                  columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
    })
    
    output$plotp_overlap <- renderPlot({
      
      
      data_selected() %>%
        ggplot(aes(p_overlap)) +
        geom_histogram(binwidth = 10, fill = '#FDE725FF', color = 'black', alpha = 0.6) +
        theme_fancy() +
        xlab('Percentage of overlapping (%)') +
        ylab('Number of genes')
      
      
      
      
    })
    
    output$plot_size <- renderPlot({
      
      req(coord_user())
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)
      
      size_cnv_query = end_coordinates - start_coordinates + 1
      
      
      if (input$select_density == 'global') {
        
      cnv_df %>%
        ggplot(aes(length_cnv, y = source)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE, size = 1.25) +
        geom_vline(aes(xintercept = size_cnv_query), linetype = 2, color = 'red', size = 1.5) +
        scale_x_log10() +
        scale_y_discrete(expand = c(0.01, 0)) +
        scale_fill_viridis_d() +
        # scale_fill_manual(values = c('#CD5C5C','#32CD32', '#32CD32')) +
        xlab('log10(CNVs size)') +
        ylab('Database') +
        theme_ridges()
      } else {
        check_cnv_df() %>%
          ggplot(aes(length_cnv, y = source)) +
          stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, show.legend = FALSE, size = 1.25) +
          geom_vline(aes(xintercept = size_cnv_query), linetype = 2, color = 'red', size = 1.5) +
          scale_x_log10() +
          scale_y_discrete(expand = c(0.01, 0)) +
          scale_fill_viridis_d() +
          xlab('log10(CNVs size)') +
          ylab('Database') +
          theme_ridges()
        
        
        
        
      }
      
      
      
    })
    
    output$tissue_gtex <- renderPlotly({
      
      req(input$input_gene_tissue)
      
      
      filtered_gene <- input$input_gene_tissue
      
      p <- gtex %>%
        filter(gene == !!filtered_gene) %>%
        ggplot(aes(reorder(tissue, -value), value)) +
        geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
        theme_fancy() +
        xlab('Tissue') +
        ylab(paste('Median TPM')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position='none') +
        ggtitle(paste('Gene expression:', filtered_gene))
      
      ggplotly(p)
      
    })
    
    output$tissue_hpa <- renderDT({
      
      req(input$input_gene_tissue)
      
      filtered_gene <- input$input_gene_tissue
      
      datatable(hpa %>% filter(gene == !!filtered_gene),
                options = list(searchHighlight = TRUE), filter = 'top', style = 'bootstrap')
      
      
      
    })
    
    
    # output$model_phenotype <- renderDT({
    #   
    #   test_tmp <- data_selected() %>% select(gene) %>% pull()
    #   mgi_test <- mgi %>% filter(gene %in% test_tmp)
    #   mgi_test <- mgi_test %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
    #     gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
    #     filter(mpo_id != '') %>%
    #     select(-delete)
    #   
    #   vector_mpo <- mgi_test %>% select(mpo_id) %>% unique() %>%
    #     mutate(description = map_chr(vector_mpo$mpo_id, function(x) termLabel(term(go, x))))
    #   
    #   mgi_test <- mgi_test %>% left_join(vector_mpo) %>% select(gene, mpo_id, description) %>% count(description)
    #   test16 <<- mgi_test
    #   
    #   
    # })
    
    
    
    output$model_genes <- renderDT({
      
      df_tmp <- model_genes_phenotype() %>% select(-entrez_id, -gene_mouse)
      
      datatable(df_tmp,
                options = list(searchHighlight = TRUE,  style = 'bootstrap'))
      
      
      
    })
    
    
    output$perc_gene <- renderPlot({
      
      validate(
        need(input$dgenes_rows_selected != '', "Please, select a gene in the datatable.")
      )
      
      red_line <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(pLI) %>% pull()
      symbol_chosen <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(gene) %>% pull()
      
      hgcn_genes %>% 
        ggplot(aes(pLI)) + 
        geom_histogram(binwidth = 0.1) +
        theme_fancy() +
        scale_fill_viridis_d() +
        geom_vline(xintercept = red_line, color = 'red', alpha = 0.6, type = 'dashed') +
        ggtitle(paste('Histogram pLI scores - Gene:', symbol_chosen))
      
    })
    
    output$dot_comparison <- renderPlotly({
      
      validate(
        need(input$dgenes_rows_selected != '', "Please, select a gene in the datatable.")
      )
      
      name_gene_filtered <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(gene) %>% pull()
      
      p <- data_selected() %>%
        
        ggplot(aes(pLI, rvis)) +
        geom_point(col = "steelblue") +
        theme_fancy() +
        gghighlight(gene == name_gene_filtered)
      
      p <- ggplotly(p)
      p <- plotly_build(p)
      
      p$x$data[1][[1]]$text <- paste0('Gene symbol: ', data_selected() %>% pull(gene), 
                                      "<br>",
                                      'pLI: ', data_selected()  %>% pull(pLI),
                                      "<br>",
                                      'RVIS: ', data_selected()  %>% pull(rvis)
                                      
      )
      
      p$x$data[[2]]$text <- paste0('Gene symbol: ', data_selected()  %>% slice(input$dgenes_rows_selected) %>% pull(gene), 
                                   "<br>",
                                   'pLI: ', data_selected() %>% slice(input$dgenes_rows_selected) %>% pull(pLI),
                                   "<br>",
                                   'RVIS: ', data_selected() %>% slice(input$dgenes_rows_selected)  %>% pull(rvis)
      )
      
      test1898 <<- p
      p
      
      
    })
    
    
    
    
    
    
    output$data <- renderTable({
      mtcars[, c("mpg", input$variable), drop = FALSE]
    }, rownames = TRUE)
    
    
    
    plot_chrom_react <- reactive({
      
      
      # req(input$start_analysis > 0)
      req(coord_user())
      
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- paste0('chr', coord_user()[3])
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)
      
      # test1916 <<-  coord_user()[1]
      # test1917 <<- coord_user()[2]
      # test1918 <-  paste0('chr', coord_user()[3])
      
      # data_input <- data_selected()
      # plotKaryotype()
      
      # ideoTrack <- IdeogramTrack(genome="hg19", chromosome= input$input_chr)
      # plotTracks(ideoTrack, from= input$int_start , to= input$int_end, showBandId=TRUE,
      #            cex.bands=0.5)
      # 
      ideoTrack <- IdeogramTrack(genome="hg19", chromosome= chrom_coordinates)
      plotTracks(ideoTrack, from= start_coordinates , to = end_coordinates, showBandId=TRUE,
                 cex.bands=0.5)
      
      # plotKaryotype(chromosomes = input_chr, plot.type = 2) %>%
      # kpDataBackground(data.panel = 1)
      # kpAddBaseNumbers() %>%
      # kpRect(chr= input_chr, x0= data_input$start_position, x1= data_input$end_position, y0=0.2, y1=0.4)
      
      
    })
    
    
    fisher_running <- reactive({
      
      req(input$start_analysis > 0)
      
      
      list_panel_names <- panel_total %>% select(Level4) %>% distinct() %>% pull()
      input_genes <- data_selected() %>% select(gene) %>% pull()
      
      fisher_result <- tibble()
      
      for (i in 1:length(list_panel_names)) {
        # print(i)
        input_test <- matrix(rep(NA, 4), ncol = 2, dimnames = list(c('in_panel', 'no_panel'), c('col1', 'col2') ))
        
        input_test[1,][2] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% pull(gene) %in% input_genes %>% sum()
        input_test[1,][1] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% nrow() - input_test[1,][2] 
        input_test[2,][2] <- length(input_genes) - input_test[1,][2] 
        input_test[2,][1] <- 19146 - input_test[2,][2] - input_test[1,][1] - input_test[1,][2]
        
        name_genes <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% pull(gene)
        name_genes <- name_genes[name_genes %in% input_genes]
        if (length(name_genes) == 0) {
          name_genes <- '-'
        }
        
        tmp_tibble <- tibble(name_panel = list_panel_names[i], 
                             p_value = fisher.test(input_test, alternative = 'greater')$p.value,
                             genes  = name_genes,
                             gene_ratio = paste(as.character(input_test[1,][2]), '/',
                                                as.character(input_test[1,][2] + input_test[1,][1]))
        )
        fisher_result <- rbind(fisher_result, tmp_tibble)
      }
      
      fisher_result %>% arrange(p_value) %>% mutate(p_value = round(p_value, 4))

    })
    
    
    output$df_fisher <- renderDataTable({
      
      
      tmp_df <- fisher_running() %>% select(-genes)
      datatable(tmp_df,
                selection = 'single',
                colnames = c('Panel name', 'p.value', 'Gene ratio'))
      
    })
    
    output$df_fisher_selected <- renderText({
      
      test88777 <<- input$df_fisher_rows_selected
      test88888 <<- fisher_running()
      validate(
        need(input$df_fisher_rows_selected != '', "Please, select a clinical panel in the datatable.")
      )
      
      name_panel_to_filter <- fisher_running() %>% slice(input$df_fisher_rows_selected) %>% pull(name_panel)
      tmp_df <- fisher_running() %>% filter(name_panel == name_panel_to_filter) %>% select(genes) %>%
        filter(!str_detect(genes, '-'))
      
      validate(
        need(nrow(tmp_df) != 0, "Not genes found.")
      )
      
      # datatable(tmp_df)
      test000 <<- tmp_df
      tmp_df %>% pull(genes) %>% paste(collapse = ', ')
      
    })
    
    fisher_running_two <- reactive({
      
      req(input$start_analysis > 0)
      
      
      name_dbs <- c('omim', 'clinvar', 'haplo', 'triplo', 'dev', 'fda', 'gwas', 'essent')
      input_genes <- data_selected() %>% select(gene) %>% pull()
      # input_genes <- test2019 %>% select(gene) %>% pull()
      fisher_result <- tibble()
      
      for (i in 1:length(name_dbs)) {
        # i <- 8
        print(i)
        input_test <- matrix(rep(NA, 4), ncol = 2, dimnames = list(c('in_panel', 'no_panel'), c('col1', 'col2') ))
        
        input_test[1,][2] <- hgcn_genes %>% filter(get(name_dbs[i]) == 'Yes') %>% pull(gene) %in% input_genes %>% sum()
        input_test[1,][1] <- hgcn_genes %>% filter(get(name_dbs[i]) == 'Yes') %>% nrow() - input_test[1,][2] 
        input_test[2,][2] <- length(input_genes) - input_test[1,][2] 
        input_test[2,][1] <- 19146 - input_test[2,][2] - input_test[1,][1] - input_test[1,][2]
        test312321 <<- input_test
        tmp_tibble <- tibble(name_panel = name_dbs[i], 
                             p_value = fisher.test(input_test, alternative = 'greater')$p.value,
                             gene_ratio = paste(as.character(input_test[1,][2]), '/',
                                                as.character(input_test[1,][2] + input_test[1,][1]))
        )
        test44444 <<- tmp_tibble
        fisher_result <- rbind(fisher_result, tmp_tibble)
      }
      
      fisher_result %>% arrange(p_value) %>% mutate(p_value = round(p_value, 4))
      
    })
    
    output$df_fisher_two <- renderDataTable({
      
      datatable(fisher_running_two(), colnames = c('Database name', 'p.value', 'Gene ratio'))
      
    })
    

    
    running_wilcoxon <- reactive({
      
      mtcars
      
      
      
      
    })
    
    output$df_wilcoxon <- renderDataTable({
      
      
      
      
      
      
    })
    
    
    output$df_fisher_selected_genes <- renderDataTable({
      
      input$input$dgenes_rows_all 
      
      datatable(fisher_running())
      
    })
    
    
    
    output$plot_chrom <- renderPlot({
      
      # req(nrow(data_selected() > 0))
      
      test1000 <<-coord_user()
      
      plot_chrom_react()
      
      
    })
    
    output$choose_geno_karyo1 <- renderUI({
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        numericInput(
          inputId = "int_start",
          label = "Genomic interval - Start",
          value = 34813719)
        
      } else {
        
        karyotype_filtered <- as.list(chromPlot::hg_cytoBandIdeo %>% filter(Chrom == coord_user()[3]) %>% select(Name))
        
        # selectizeInput(inputId = 'input_karyotype', label = 'Karyotype', choices = karyotype_filtered,
        #                selected = NULL, multiple = FALSE,
        #                options = NULL)
        
        pickerInput(
          inputId = "input_karyotype",
          label = "Karyotype", 
          choices = karyotype_filtered,
          options = list(
            `live-search` = TRUE)
        )
        
      }
      
    })
    
    
    
    output$n_genes <- renderUI({
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        name_region <- paste0('chr',coord_user()[3], ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(coord_user()[3], input$input_karyotype)
        
      }

      test1111111 <- data_selected()
      if (is.null(input$dgenes_rows_all)) {
        final_number <- test1111111
      } else {
        final_number <- test1111111[input$dgenes_rows_all,]
      }
      test1111111 <- test1111111[input$dgenes_rows_all,]
      test787 <<- input$dgenes_rows_all
      test2121 <<- test1111111
      
      tablerStatCard(
        value = nrow(final_number),
        title = HTML(paste0("Number of genes<br/>", name_region)),
        # trend = 19192,
        width = 12
      )

    })
    
    output$score_rf <- renderUI({
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)
      
      length_input_cnv <- end_coordinates - start_coordinates + 1
      tmp_df <- data_selected()
      
      input_model <- tibble(n_genes = NA, pli = NA, rvis = NA, omim = NA, 
                            length_cnv =  length_input_cnv)
      
      if (nrow(tmp_df) == 0) {
        
        input_model$n_genes[1]  <- 0
        input_model$pli[1]  <- 0
        input_model$rvis[1]  <- 0
        input_model$omim[1] <- 0
        
      } else {
        
        # Variable: number of genes disrupted
        input_model$n_genes[1]  <- nrow(tmp_df)
        # Variable: sum pLI score
        input_model$pli[1]  <- tmp_df %>% pull(pLI) %>% sum(na.rm = TRUE)
        # Variable: sum rvis score
        input_model$rvis[1]  <- tmp_df %>% pull(rvis) %>% sum(na.rm = TRUE)
        # Variable: nº genes included in OMIM
        tmp_omim<- tmp_df %>% count(omim) %>% filter(omim == 'Yes') %>% pull(n)
        if (length(tmp_omim) == 0) {
          input_model$omim[1] <- 0
        } else {
          input_model$omim[1] <- tmp_omim
          
        }
      }
      
      test9866 <<- input_model
  
      score_predicted <- predict(model1, input_model, type = "prob") %>% 
        as_tibble() %>% pull(decipher) %>% as.numeric() * 100
      
      tablerStatCard(
        value = score_predicted,
        title = 'Score pathogenicty (0-100) (JUST AN EXAMPLE)',
        # trend = 19192,
        width = 12
      )
      
    })
    
    # output$n_snv <- renderUI({
    #   
    #   tablerStatCard(
    #     value = nrow(reading_snv_file()),
    #     title = "Number of SNVs",
    #     # trend = 19192,
    #     width = 12
    #   )
    #   
    # })
    
    
    
    
    
    output$ref_user_filter_genes <- renderUI({
      
     req(length(input$dgenes_rows_all) != nrow(data_selected()))

      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        name_region <- paste0('chr',coord_user()[3], ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(coord_user()[3], input$input_karyotype)
        
      }
      
      tablerInfoCard(
        width = 12,
        value = paste0(length(input$dgenes_rows_all), '/', nrow(data_selected()), " genes"),
        status = "warning",
        icon = "crop",
        description =  'Filtered genes'
      )
      
      
    })
    
    output$ref_user_filter_genes2 <- renderUI({
      
      req(length(input$dgenes_rows_all) != nrow(data_selected()))
      
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        name_region <- paste0('chr',coord_user()[3], ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(coord_user()[3], input$input_karyotype)
        
      }
      
      tablerInfoCard(
        width = 12,
        value = paste0(length(input$dgenes_rows_all), '/', nrow(data_selected()), " genes"),
        status = "warning",
        icon = "crop",
        description =  'Filtered genes'
      )
      
      
    })
    
    
    output$genes_enhancers_selected <- renderUI({
      
      # req(length(input$dgenes_rows_all) != nrow(data_selected()))
      
      # if (input$input_geno_karyo == 'Genomic coordinates') {
      #   
      #   name_region <- paste0('chr',input$input_chrom, ':', input$int_start, '-', input$int_end)
      #   
      # } else {
      #   name_region <- paste0(input$input_chrom, input$input_karyotype)
      #   
      # }
      # 
      
      da
      
      tablerInfoCard(
        width = 12,
        value = paste0(length(input$dgenes_rows_all), '/', nrow(data_selected()), " genes"),
        status = "warning",
        icon = "crop",
        description =  'Filtered genes'
      )
      
      
    })
    
    output$ref_user_genes <- renderUI({
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        name_region <- paste0('chr', coord_user()[3], ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(coord_user()[3], input$input_karyotype)
        
      }
      
      tablerInfoCard(
        width = 12,
        value = paste0(nrow(data_selected()), " genes"),
        status = "success",
        icon = "database",
        description =  name_region
      )
      
      
    })
    
    output$ref_user_length <- renderUI({
      
       req(coord_user())
      # 
      # 
      # if (input$input_geno_karyo == 'Genomic coordinates') {
      #   
      #   name_region <- paste0('chr',input$input_chrom, ':', input$int_start, '-', input$int_end)
      #   
      # } else {
      #   name_region <- paste0(input$input_chrom, input$input_karyotype)
      #   
      # }
      # 
      # 
      # if (input$input_geno_karyo == 'Genomic coordinates') {
      #   
      #   length_region <-  round((input$int_end - input$int_start + 1) / 1e6, 2)
      #   
      # } else {
      #   
      #   # name_region <- paste0(input$input_chrom, input$input_karyotype)
      #   length_region <-  round((input$int_end - input$int_start + 1) / 1e6, 2)
      # }
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)
      
      
      test3133 <<- coord_user()[1]
      test3134 <<- coord_user()[2]
      
      
      length_region <- end_coordinates - start_coordinates + 1
      length_region <- round(length_region / 1e6, 2)
      
      tablerInfoCard(
        width = 12,
        value =  paste(length_region, 'Mb'),
        status = "primary",
        icon = "database",
        description =  'Length of the genomic region'

      )
      
      
    })
    
    
    output$n_genes_added <- renderUI({
      
      req(nrow(data_selected_enhancers()) > 0)
      
    
      tablerInfoCard(
        width = 12,
        value =  paste0('+', nrow(data_selected_enhancers()), ' genes'),
        status = "info",
        icon = "database",
        description =  'Target-genes enhancers'
        
      )
      
      
    })
    
    
    lncrna_raw <- reactive({
      
      req(input$start_analysis > 0)
      
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      data_tmp <- lncrna_coord %>% filter(chrom == chrom_coordinates) %>%
        mutate(keep = 0)
      
      for (i in 1:nrow(data_tmp)) {
        data_tmp$keep[i] <- c(data_tmp$start[i], data_tmp$end[i]) %overlaps% c(start_coordinates, end_coordinates)
      }
      
      data_tmp <- data_tmp %>% filter(keep == 1) %>% select(id) %>% distinct() %>% pull(id)
      
      
    })
    
    output$lncrna_df <- renderDataTable({
      
      lncrna_selected <- lncrna_raw()
      
      tmp_lncrna <- lncrna %>% filter(id %in% lncrna_selected) %>% select(-genomic_class)
      
      datatable(tmp_lncrna, rownames = FALSE, 
                colnames = c('ID', 'Chromosome', 'Dynamic', 'Ensembl75 ID', 'Name', 'Conservation', 'Nearest coding gene'),
                options = list(
                  pageLength = 5, autoWidth = TRUE, list(searchHighlight = TRUE)))
      
      
      
    })
    
    
    output$n_lncrna <- renderUI({
      
      tablerStatCard(
        value =  length(lncrna_raw()),
        title = "Number of lncRNA disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    prev_enhancer <- reactive({
      
      req(input$start_analysis > 0)
      
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      data_tmp <- df_enhancers %>% filter(chrom == chrom_coordinates) %>%
        mutate(keep = 0)
      
      for (i in 1:nrow(data_tmp)) {
        data_tmp$keep[i] <- c(data_tmp$start[i], data_tmp$end[i]) %overlaps% c(start_coordinates, end_coordinates)
      }
      
      data_tmp <- data_tmp %>% 
        # mutate(keep = c(start, end) %overlaps% c(start_coordinates, end_coordinates)) %>%
        filter(keep == 1) %>% 
        select(-keep)
      
      data_tmp <- data_tmp %>% 
        left_join(hgcn_genes %>% select(gene, chrom, start_position, end_position, chrom) %>% 
                     rename(chrom_gene = chrom, start_gene = start_position, end_gene = end_position)) %>%
        rowwise() %>%
        mutate(inside_cnv = c(start_gene, end_gene) %overlaps% c(start_coordinates, end_coordinates)) %>%
        mutate(inside_cnv = if_else(chrom_gene == chrom_coordinates, inside_cnv, FALSE)) %>%
        select(-chrom_gene, -start_gene, -end_gene) %>%
        distinct()
        
        test1946 <<- data_tmp
        
      # get_coord_genes <- hgcn_genes %>% 
      #   filter(gene %in% (data_tmp %>% pull(gene))) %>%
      #   mutate(inside_cnv = c(start, end) %overlaps% c(start_coordinates, end_coordinates)) %>%
      #   select(gene, inside_cnv)
      # 
      #   test1946 <<- get_coord_genes
      #   
      #   data_tmp <- data_tmp %>% left_join(get_coord_genes)
        

      data_tmp
      
      
    })
    
    # output$n_enhancer_inside <- renderUI({
    #   
    #   test1936 <<- prev_enhancer()
    #   
    #   list_genes_total <- prev_enhancer() %>% select(gene) %>% distinct() %>% pull(gene)
    #   data_tmp <- data_selected() %>% filter(gene %in% list_genes_total) %>% pull(gene)
    #   
    #   tablerStatCard(
    #     value =  paste(length(data_tmp), length(list_genes_total), sep = '/'),
    #     title = "Number of target genes inside CNV",
    #     width = 12
    #   )
    # })
    
    
    
    output$n_enhancer <- renderUI({
      
      data_tmp <- prev_enhancer() %>% select(id) %>% distinct() %>% pull(id)
      
      
      tablerStatCard(
        value =  length(data_tmp),
        title = "Number of Enhancers disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    
    redundancy_enhancers <- reactive({
      
      data_tmp <- prev_enhancer() %>% select(gene) %>% distinct() %>% pull(gene)
      df_output <- df_enhancers %>% filter(gene %in% data_tmp) %>% count(gene)
      test25 <<- df_output
    })
    
    # output$redund_n_enhancer <- renderUI({
    #   
    #   test24 <<-  redundancy_enhancers() 
    #   
    #   data_tmp <- redundancy_enhancers() %>% filter(n == 1)
    # 
    #   tablerStatCard(
    #     value =  nrow(data_tmp),
    #     title = "Number of genes whose have one enhancer and it is disrupted",
    #     # trend = -10,
    #     width = 12
    #   )
    # 
    # })
    
    output$histo_enhancer <- renderPlot({
      
      validate(
        need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
      )
      
      score_filtered <- prev_enhancer() %>% slice(input$df_enhancer_rows_selected) %>% select(phast100) %>% pull()

      prev_enhancer() %>% ggplot(aes(phast100)) + 
        geom_histogram() + 
        theme_fancy() +
        geom_vline(xintercept = score_filtered, color = 'red')



    })
    
    
    output$switch_enhancers <- renderUI({
      
      
      validate(
        need(nrow(prev_enhancer()) > 0, "Need a region with at least one enhancer.")
      )
      
      
      switchInput(
        inputId = "enhancers_on_off",
        label = "Add target genes to analysis?",
        inline = TRUE,
        width = 'auto',
        # status = "warning",
        value = FALSE
        # right = TRUE
      )
      
    })
    
    output$df_enhancer <- renderDataTable({
      
      
      df_tmp <- prev_enhancer() 
      df_tmp <- df_tmp %>% select(-score_enh, -score_gene)
    
      
      datatable(df_tmp, rownames = FALSE, filter = 'top', selection = 'single',
                colnames = c('Gene', 'ID enhancer', 'Chrom',
                              'Start',  'End', 'Phast100way', 'Phast46way Placental', 'Phast46way Primates','o/e gnomad', 
                             'Gene mapped in CNV'),
                options = list(
                  pageLength = 5, 
                  # autoWidth = TRUE,
                  # style = 'bootstrap',
                  list(searchHighlight = TRUE)
                  # selection = 'single'
                  # columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
      
      
    })
    
    output$n_tads <- renderUI({
      
      req(input$start_analysis > 0)
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      
      n_tads <- check_tads(chrom_coordinates, start_coordinates, end_coordinates, tad )
      n_tads <- nrow(n_tads)
      test931323 <<- n_tads
      
      if (is.null(n_tads)) {
        n_tads <- 0
      } else {
        n_tads <- as.double(n_tads)
      }
      
      tablerStatCard(
        value =  n_tads,
        title = "Number of TADs disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    
output$n_filtered_enhancers <- renderUI({
      
  
  if (is.null(input$df_enhancer_rows_all)) {
    n_enhancers <- prev_enhancer()
  } else {
    n_enhancers <- prev_enhancer()[input$df_enhancer_rows_all,]
  }
  
  n_total_enhancers <- prev_enhancer() %>% nrow()
  
  n_enhancers <- n_enhancers %>% nrow()
    
    tablerInfoCard(
      width = 12,
      value =  paste0(n_enhancers, '/', n_total_enhancers, ' target genes'),
      status = "warning",
      icon = "database",
      description =  'Filtered enhancers'
      
    )
  
})
    
    output$df_tads <- renderDataTable({
      
      req(input$start_analysis > 0)
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      
      n_tads <- check_tads(chrom_coordinates, start_coordinates, end_coordinates, tad )
      
      validate(
        need(!is.null(nrow(n_tads)), "0 TADs found.")
      )
      
      datatable(n_tads,
                colnames = c('ID', 'Chromosome', 'Start', 'End'))
    })
    
    output$n_dev <- renderUI({
      
      
      n_dev_yes <- data_selected() %>% filter(dev == 'Yes') %>% nrow()
      n_total <- nrow(data_selected())
      tablerStatCard(
        value =  paste(n_dev_yes, n_total, sep = '/'),
        title = "Developmental disorder genes",
        # trend = -10,
        width = 12
      )
    })
    
    output$n_clinvar <- renderUI({
      
      
      n_clinvar_yes <- data_selected() %>% filter(clinvar == "Yes") %>% nrow()
      n_total <- nrow(data_selected())  
      
      tablerStatCard(
        value =  paste(n_clinvar_yes, n_total, sep = '/'),
        title =  'Genes found in ClinVar',
        # trend = -10,
        width = 12
      )
    })
    
    output$n_gwas <- renderUI({
      
      
      n_gwas_yes <- data_selected() %>% filter(gwas == 'Yes') %>% nrow()
      n_total <- nrow(data_selected())
      
      tablerStatCard(
        value =  paste(n_gwas_yes, n_total, sep = '/'),
        title =  'Genes found in GWAS Catalog',
        # trend = -10,
        width = 12
      )
    })
    
    
    output$n_hi<- renderUI({
      
      
      n_hi_yes <- data_selected() %>% filter(hi <= 10) %>% nrow()
      n_total <- nrow(data_selected())  
      
      tablerStatCard(
        value =  paste(n_hi_yes, n_total, sep = '/'),
        title =  paste("Likely to exhibit haploinsufficiency", '(hi <= 10)'),
        # trend = -10,
        width = 12
      )
    })
    
    output$n_pli <- renderUI({
      
      test15 <- data_selected() %>% select(entrez_id) %>% pull()
      
      
      n_pli_yes <- data_selected() %>% filter(pLI >= 0.9) %>% nrow()
      n_total <- nrow(data_selected())  
      
      tablerStatCard(
        value =  paste(n_pli_yes, n_total, sep = '/'),
        title =  paste("Intolerant to LoF mutations", '(pLI >= 0.9)'),
        # trend = -10,
        width = 12
      )
    })
    
    output$n_omim <- renderUI({
      
      gene_id <- data_selected() %>% filter(omim == 'Yes') %>% pull(gene)
      n_total <- data_selected() %>% nrow()
      
      # morbidmap
      
      tablerStatCard(
        value =  paste(length(gene_id), n_total, sep = '/'),
        title = "OMIM genes",
        # trend = -10,
        width =  12
      )
    })
    
    output$p_pli <- renderUI({
      
      tablerStatCard(
        value =  hgcn_genes %>% mutate(p_pli = ntile(pLI, 100)) %>% filter(gene == 'AADACL4') %>% select(p_pli) %>% pull(),
        title = "Percentile pLI score",
        # trend = -10,
        width = 12
      )
    })
    
    output$p_rvis <- renderUI({
      
      tablerStatCard(
        value =  hgcn_genes %>% mutate(p_rvis = ntile(rvis, 100)) %>% filter(gene == 'AADACL4') %>% select(p_rvis) %>% pull(),
        title = "Percentile RVIS score",
        # trend = -10,
        width = 12
      )
    })
    
    
    
    output$n_genes_pli <- renderUI({
      
      value_input <- data_selected() %>% filter(pLI >= 0.9) %>% nrow()
      tablerInfoCard(
        value = value_input,
        status = "danger",
        icon = "dollar-sign",
        description = ">0.9 pLI score",
        width = 12
      )
      
      
    })
    
    
    
    output$choose_geno_karyo2 <- renderUI({
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        numericInput(
          inputId = "int_end",
          label = "Genomic interval - End",
          value = 36278623)
      }
    })
    
    model_genes_phenotype <- reactive({
      
      
      go <- Ontology("hp")
      
      test_tmp <- data_selected() %>% select(gene) %>% pull()
      mgi_tmp <- mgi %>% filter(gene %in% test_tmp)
      mgi_tmp <- mgi_tmp %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
        gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
        filter(mpo_id != '') %>%
        select(-delete)
      
      vector_mpo <- mgi_tmp %>% select(mpo_id) %>% unique() %>%
        mutate(description = map_chr(mpo_id, function(x) termLabel(term(go, x)))) %>%
        mutate(description = str_remove(description, ' phenotype'))
      
      mgi_tmp <- mgi_tmp %>% left_join(vector_mpo)
      
      mgi_tmp
      
      
      
    })
    
    
    output$agg_model <- renderEcharts4r({
      
      model_genes_phenotype() %>% count(description) %>% arrange(n) %>%
        e_charts() %>% 
        e_treemap(description, description, n) %>%
        e_tooltip(trigger = "axis") %>%
        e_title("")
      
      
    })
    
    output$funnel_genes <- renderEcharts4r({
      
      genes_selected <- paste(nrow(data_selected()), 'genes found in CNV') 
      
      test1000 <<- input$dgenes_rows_all
      
      if (nrow(data_selected()) == length(input$dgenes_rows_all)) {
        genes_filtered <-data_selected()
        funnel <- data.frame(stage = c("19,146 genes", genes_selected), value = c(1, 0.5))
        
      } else {
        genes_filtered <-  paste(nrow(data_selected()[input$dgenes_rows_all,]), 'genes filtered')
        funnel <- data.frame(stage = c("19,146 genes", genes_selected, genes_filtered), value = c(1, 0.5, 0.25))
      }

      funnel %>% 
        e_charts() %>% 
        e_funnel(value, stage) %>% 
        e_title("")
      
    })
    
    
    # output$info <- renderUI({
    #   tablerInfoCard(
    #     width = 12,
    #     value = paste0(input$totalStorage, "GB"),
    #     status = "success",
    #     icon = "database",
    #     description = "Total Storage Capacity"
    #   )
    # })
    
    
    # output$heatmap <- renderPlotly({
    #   
    #   # a <- matrix(NA, nrow = 2, ncol = 100)
    #   # a[1,] <- hgcn_genes$pLI[1:100]
    #   # a[2,] <- hgcn_genes$rvis[1:100]
    #   # a[1,][as.numeric(which(is.na(a[1,])))] <- 0
    #   # a[2,][as.numeric(which(is.na(a[2,])))] <- 0
    #   # colnames(a) <- hgcn_genes$gene[1:100]
    #   # a <- pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE)
    #   # a
    #   n_genes <- nrow(data_selected())
    #   data_raw <- data_selected() %>% mutate(p_li = ntile(pLI, 100), p_rvis = ntile(rvis, 100),
    #                                          p_ncrvis = ntile(ncrvis, 100), p_ncgerp = ntile(ncgerp, 100))
    #   m <- matrix(NA, nrow = 4, ncol = n_genes)
    #   m[1,] <- data_raw$p_li[1:n_genes]
    #   m[2,] <- data_raw$p_rvis[1:n_genes]
    #   m[3,] <- data_raw$p_ncrvis[1:n_genes]
    #   m[4,] <- data_raw$p_ncgerp[1:n_genes]
    #   # m <- m[colSums(!is.na(m)) > 0]
    #   
    #   # heatmaply(as.matrix(m))
    #   
    #   plot_ly(
    #     x = data_selected() %>% select(gene) %>% pull(), y = c("RVIS", "pLI", 'ncRVIS', 'ncGERP'),
    #     z = m, 
    #     type = "heatmap"
    #     # width = 1200,
    #     # height = 500
    #   )
    # })
    
    
    running_enrich_go <- reactive({
      
      req(isTRUE(input$enable_func_analysis))
      
        if (input$choose_go == 'biological process') {
          go_chosen <- 'BP'
        } else if (input$choose_go == 'molecular function') {
          go_chosen <- 'MF'
        } else {
          go_chosen <- 'CC'
        }
      
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      
        
        filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
        univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character()
        
        validate(
          need(length(filtered_genes) != 0, "0 enriched terms found.")
        )
        
        go_analysis <- enrichGO(gene  = filtered_genes,
                                universe      = univ,
                                OrgDb         = org.Hs.eg.db,
                                ont           = go_chosen,
                                # pAdjustMethod = "BH",
                                pvalueCutoff  = as.numeric(input$sign_vline),
                                readable = TRUE)
        
                                # qvalueCutoff  = 0.05)
        

        go_analysis <- go_analysis %>% as_tibble()
      
        validate(
          need(nrow(go_analysis) != 0, "0 enriched terms found.")
        )
      
        go_analysis
      })
    
    
    output$func_analysis <- renderPlot({
      
      test455 <<- running_enrich_go() 
    
       running_enrich_go() %>%
          mutate(p.adjust = -log10(p.adjust)) %>%
          # arrange(desc(Count)) %>%
          slice(1:20) %>%
          ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
          geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
          scale_fill_viridis_d() +
          coord_flip() +
          xlab('') +
          ylab('-log10(p-adjusted)') +
          ggtitle('GO analysis') +
          geom_hline(yintercept =-log10(as.numeric(input$sign_vline)), color = 'red', alpha = 0.6, linetype = 'dashed', size = 2) +
          geom_text(
            aes(label = GeneRatio, y = p.adjust + 0.05),
            position = position_stack(vjust = 0.5),
            vjust = 0
          ) +
         theme_fancy() +
          theme(axis.text=element_text(size=12),
                axis.title=element_text(size=14,face="bold"))
        
        

    })
    
    
    
    output$df_enrichgo <- renderDT({
      
      df <- running_enrich_go()  %>%
        as_tibble() %>%
        filter(Count != 0) %>%
        arrange(desc(Count)) %>%
        separate(geneID, sep = '/', into = as.character(1:1000)) %>%
        gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
        select(-delete) %>%
        na.omit() %>%
        distinct() %>%
        mutate(pvalue = round(pvalue, 3)) %>%
        mutate(p.adjust = round(p.adjust, 3)) %>%
        mutate(qvalue = round(qvalue, 3))
        
      
      datatable(df)
      
      
    })
    
    running_go <- reactive({
      
      
      req(isTRUE(input$enable_group_go))
      
      
      if (input$choose_group_go == 'biological process') {
        go_chosen <- 'BP'
      } else if (input$choose_group_go == 'molecular function') {
        go_chosen <- 'MF'
      } else {
        go_chosen <- 'CC'
      }
      
      n_level <- as.numeric(input$user_level)
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      
      filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
      
      tryCatch(
        
        ggo <- groupGO(gene     = filtered_genes,
                       OrgDb    = org.Hs.eg.db,
                       ont      = go_chosen,
                       level    = n_level,
                       readable = TRUE),
        
        
        error= function(e) stop("Please, reduce the level assigned"))
      
      ggo <- ggo %>% as_tibble()
      
      
      validate(
        need(ggo %>% as_tibble() %>% select(Count) %>% sum() != 0, "0 terms found. Please reduce the level assigned")
      )
      
      ggo
      
    })
    
    
    output$plot_enrichgo <- renderUI({
      
      if (input$table_plot_go == 'Table') {
        
        DTOutput('df_enrichgo')
        
      } else {
        plotOutput('func_analysis')
        
      }
      
    })
    
    
    output$plot_df <- renderUI({
      
      if (input$user_df_plot == 'Table') {
        
        DTOutput('df_go')
        
      } else {
        plotOutput('group_go')
        
      }
      
    })
    
    
    output$group_go  <- renderPlot({
      
      test900 <<- running_go()
      test344 <<- running_enrich_go()
      running_go() %>%
        arrange(desc(Count)) %>%
        filter(Count != 0) %>%
        slice(1:10) %>%
        as_tibble() %>% 
        # mutate(p.adjust = -log10(p.adjust)) %>%
        ggplot(aes(reorder(Description, Count), Count)) +
        geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
        scale_fill_viridis_d() +
        coord_flip() +
        xlab('') +
        ylab('Number of genes') +
        geom_text(
          aes(label = GeneRatio, y = Count + 0.05),
          position = position_stack(vjust = 0.5),
          vjust = 0
        ) +
        theme_fancy() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
      
      
    })
    
    
    output$df_go  <- renderDT({
      
      
      
      test222 <<- running_go()
      
      df <- running_go() %>%
        as_tibble() %>%
        filter(Count != 0) %>%
        arrange(desc(Count)) %>%
        separate(geneID, sep = '/', into = as.character(1:1000)) %>%
        gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio) %>%
        select(-delete) %>%
        na.omit() %>%
        distinct()
      
      datatable(df)
      
    })
    
    running_path <- reactive({
      
      req(isTRUE(input$enable_path_analysis))
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
      

      validate(
        need(length(filtered_genes) != 0, "0 enriched terms found.")
      )
      
      if (input$kegg_reactome == 'Reactome') {
        
        pathway_analysis <- enrichPathway(gene= filtered_genes ,
                                          pvalueCutoff = as.numeric(input$sign_vline_path), 
                                          universe = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
                                          readable= TRUE)  %>% 
          as_tibble()
      } else {
        pathway_analysis <- clusterProfiler::enrichKEGG(gene= filtered_genes ,
                                                        pvalueCutoff= as.numeric(input$sign_vline_path),
                                                        universe = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character()) %>% 
          as_tibble()
      }
      
      
      validate(
        need(nrow(pathway_analysis) != 0, "0 enriched terms found.")
      )      
      
      pathway_analysis %>%
        mutate(pvalue = round(pvalue, 3)) %>%
        mutate(p.adjust = round(p.adjust, 3)) %>%
        mutate(qvalue = round(qvalue, 3))
      
    })
    
    
    
    output$ui_path <- renderUI({
      
      if (input$table_path_go == 'Table') {
        
        DTOutput('df_path')
        
      } else {
        plotOutput('func_pathways')
        
      }
      
    })
    
    
    # running_do <- reactive({
    #   
    #   
    #   
    #   
    #   
    # })
    
    
    

    
    output$df_path  <- renderDT({
      
    test9991 <<- running_path()
      df <- test9991  %>%
        # filter(Count != 0) %>%
        arrange(desc(Count)) %>%
        separate(geneID, sep = '/', into = as.character(1:1000)) %>%
        gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
        select(-delete) %>%
        # na.omit() %>%
        distinct()
      
      datatable(df)

    })
    
    
    
      
    
    
    output$func_pathways  <- renderPlot({
      
      test800 <<- running_path()
      
      running_path() %>%
        mutate(p.adjust = -log10(p.adjust)) %>%
        ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
        geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
        scale_fill_viridis_d() +
        coord_flip() +
        xlab('') +
        ylab('-log10(p-adjusted)') +
        ggtitle('Enriched pathways') +
        geom_hline(yintercept =-log10(as.numeric(input$sign_vline_path)), color = 'red', alpha = 0.6, linetype = 'dashed', size = 2) +
        geom_text(
          aes(label = GeneRatio, y = p.adjust + 0.05),
          position = position_stack(vjust = 0.5),
          vjust = 0
        ) +
        theme_fancy() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
    })
    
    output$ui_do <- renderUI({
      
      if (input$table_plot_do == 'Table') {
        
        DTOutput('df_do')
        
      } else {
        
        plotOutput('func_do')
        
      }
      
    })
    
    running_do  <- reactive({
      
      # enrichNCG
      # enrichDGN
      # gseNCG
      # gseDGN
      
      req(isTRUE(input$enable_do_analysis))
      
      
      if (is.null(input$dgenes_rows_all)) {
        df_genes <- data_selected()
      } else {
        df_genes <- data_selected()[input$dgenes_rows_all,]
      }
      
      filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
      validate(
        need(length(filtered_genes) != 0, "0 enriched terms found.")
      )

      enrich_dgn <- enrichDO(gene  = filtered_genes,
                                                    universe      = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
                                                    # pAdjustMethod = "BH",
                                                    pvalueCutoff  = as.numeric(input$pvalue_do),
                                                    readable = TRUE)
      
      test444 <<- enrich_dgn
      
      validate(
        need(nrow(enrich_dgn) != 0, "0 enriched terms found.")
      )
      
      enrich_dgn

    })
    
    
    output$df_do <- renderDataTable({
      
      test883 <<- running_do()  
      
      df <- running_do()  %>%
        as_tibble() %>%
        filter(Count != 0) %>%
        arrange(desc(Count)) %>%
        separate(geneID, sep = '/', into = as.character(1:1000)) %>%
        gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
        select(-delete) %>%
        na.omit() %>%
        distinct() %>%
        mutate(pvalue = round(pvalue, 3)) %>%
        mutate(p.adjust = round(p.adjust, 3)) %>%
        mutate(qvalue = round(qvalue, 3))

      datatable(df)
      
    })

output$func_do  <- renderPlot({
  
  test3222 <<- running_do()
  
    cnetplot(running_do(), foldChange= hgcn_genes$entrez_id, readable = TRUE)
})
    
    output$table_all_genes  <-  renderReactable({
      
     reactable(hgcn_genes,
               filterable = TRUE,
               selection = "single", 
               selectionId = "selected",
               
               columns = list(chrom = colDef(name = "Chromosome"),
                              start_position = colDef(name = "Start"),
                              end_position = colDef(name = "End"),
                              location = colDef(name = "Location"),
                              gene = colDef(name = "Gene"),
                              haplo = colDef(name = "Haploinsufficiency",
                                             cell = function(value) {
                                               
                                               if (value == 0) "\u2718" else "\u2713"
                                             }),
                              triplo = colDef(name = "Triplosensitivity",
                                              cell = function(value) {
                                                
                                                if (value == 0) "\u2718" else "\u2713"
                                              }),
                              dev = colDef(name = "Developmental disorder gene",
                                           cell = function(value) {
                                             
                                             if (value == 0) "\u2718" else "\u2713"
                                           }),
                              fda = colDef(name = "FDA-approved drug targets",
                                           cell = function(value) {
                                             
                                             if (value == 0) "\u2718" else "\u2713"
                                           }),
                              clinvar = colDef(name = "ClinVar genes",
                                               cell = function(value) {
                                                 
                                                 if (value == 0) "\u2718" else "\u2713"
                                               }),
                              ccr = colDef(name = "CCR"),
                              ncrvis = colDef(name = "Ncrvis"),
                              ncgerp = colDef(name = "Ncgerp"),
                              hi = colDef(name = "HI index", 
                                          format = colFormat(suffix = " %", digits = 1))
                              # cell = function(value) {
                              #   if (is.na(value)) {
                              #     classes <- "tag num-low"
                              #   } else if (value >= 0.9) {
                              #     classes <- 'tag num-high'
                              #   } else  {
                              #     classes <- "tag num-high"
                              #   }
                              #   value <- format(value, nsmall = 1)
                              #   span(class = classes, value)
                              # })
               ))
      
    })
    
    
    output$button_download <- downloadHandler(
      
      filename = 'report.html',
      
      content = function(file) {
        
        params <- list(start = input$int_start,
                       end = input$int_end,
                       chrom = input$input_chrom,
                       name = input$name_report,
                       age = input$age_report,
                       sex = input$sex_report,
                       comment = input$comment_report)
                       # gene_content = data_selected())
        
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        rmarkdown::render(tempReport, output_file = file)
                          params = params
                          # envir = new.env(parent = globalenv())
      })
    
    
    # reading_snv_file <- reactive({
    #   
    #   # conditional of errors!
    #   # req(input$snv_yes_no == 'Yes')
    #   req(input$file_snv)
    #   
    #   file1 <- input$file_snv
    #   data1 <- read.table(file1$datapath, sep = '\t', header = FALSE) %>%
    #             rename(chrom = V1, start = V2, end = V3)
    # })
    
    output$test_reading_snv_file <- renderDT({
      
      datatable(reading_snv_file())
      
    })
    
    df_overlap_cnvs_running <- reactive({
      
      req(input$start_analysis > 0)
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      
      tmp_df <-  get_perc_overlap(check_cnv_df() %>% 
                                    rename(start_position = start, end_position = end), start_coordinates, 
                                  end_coordinates)

      tmp_df
    })
    
    intersection_running <- reactive({

      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)

      df_pathogenic <- df_overlap_cnvs_running() %>% filter(source == 'decipher')

      cnv_input <- GRanges(
        seqnames= chrom_coordinates,
        ranges= IRanges(start= start_coordinates, end = end_coordinates)
      )

      cnvs_pathogenic <- GRanges(
        seqnames= df_pathogenic$chrom,
        ranges= IRanges(start= df_pathogenic$start_position, end = df_pathogenic$end_position)
      )

      gr_intersect <- intersect(cnv_input, cnvs_pathogenic) %>% as_tibble()
      
      gr_intersect <- gr_intersect %>%
        mutate(id = as.factor(seq(1:n())))

      test711 <<- gr_intersect

      gr_intersect


    })

    output$df_intersection <- renderDataTable({
      
      df_tmp <- intersection_running() %>%
        select(id, start, end)

      datatable(df_tmp, rownames = FALSE, selection = 'single',
                options = list(dom = 't'))

    })

    output$plot_intersection <- renderPlot({

      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      
      start_coordinates <- as.numeric(start_coordinates)
      end_coordinates <- as.numeric(end_coordinates)

      
      df_tmp <-intersection_running()

      df_genes <- data_selected() %>% mutate(keep = NA) %>%
        rowwise() %>%
        mutate(keep = c(start_position, end_position) %overlaps% c(df_tmp$start, df_tmp$end)) %>%
        filter(keep == TRUE) %>%
        select(-keep) %>%
        ungroup() %>%
        mutate(pos = round((end_position + start_position)/2), 0) %>%
        mutate(id = rep('gene(s)', n())) %>%
        select(id, pos)

      df_total <- intersection_running() %>%
        select(id, start, end) %>%
        gather('rm', 'pos', -id) %>%
        select(-rm) %>%
        rbind(df_genes) %>%
        mutate(id_color = as.factor(if_else(id == 'gene(s)', 'steelblue', 'red')))
      
        test1453 <<- df_total

      ggplot() +
        geom_point(data = df_total, aes(pos, id, fill = id_color), color = 'black', shape = 21) +
        geom_path(data =df_total %>% filter(id != 'gene(s)'), aes(pos, id, group = id, color = id_color)) +
        coord_cartesian(xlim = c(start_coordinates, end_coordinates )) +
        ylab('Id intersections') +
        xlab(paste0('chr', coord_user()[3], ':', coord_user()[1], '-', coord_user()[2])) +
        theme_fancy()

    })

    
    df_overlap_cnvs_running_download_patho <- reactive({
      

      tmp_df <-  df_overlap_cnvs_running() %>% filter(source == 'decipher')
      
      tmp_df
    })
    
    df_overlap_cnvs_running_download_nonpatho <- reactive({
      

      tmp_df <-  df_overlap_cnvs_running() %>% filter(source != 'decipher')
      
      tmp_df
    })
    
    output$plot_cnv_bar <- renderPlot({
      
      
      start_pos <- as.numeric(coord_user()[1])
      end_pos <- as.numeric(coord_user()[2])  
      

      
      df_nonpathogenic <- df_overlap_cnvs_running() %>% filter(source != 'decipher')
      df_pathogenic <- df_overlap_cnvs_running() %>% filter(source == 'decipher')
      
      test8231 <<- start_pos
      test8999 <<- end_pos
      test88881 <<- df_nonpathogenic
      test88882 <<- df_pathogenic
      
      # start_pos <- test8231
      # end_pos <- test8999
      # df_nonpathogenic <- test88881
      # df_pathogenic <- test88882
      
      
      interval_values <- seq(from = start_pos, to = end_pos, length.out = 200)
      
      df_interval <- matrix(interval_values, ncol = 2, byrow = TRUE)
      colnames(df_interval) <- c('start', 'end')
      df_interval <- as_tibble(df_interval)
      
      query <- IRanges(df_interval$start, df_interval$end)
      
      # Evaluate pathogenic CNVs
      subject_patho <- IRanges(df_pathogenic$start_position, df_pathogenic$end_position)
      hits_intervals_patho <- countOverlaps(query, subject_patho)
      
      # Evaluate pathogenic CNVs
      subject_nonpatho<- IRanges(df_nonpathogenic$start_position, df_nonpathogenic$end_position)
      hits_intervals_nonpatho <- countOverlaps(query, subject_nonpatho)
      
      
      df_interval <- df_interval %>% mutate(n_patho = hits_intervals_patho, 
                                            n_nonpatho = hits_intervals_nonpatho) %>%
        mutate(id = row_number()) %>%
        gather('category', 'n_overlap', -start, -end, -id)
      
      
      df_interval %>%
        mutate(category = if_else(category == 'n_patho', 'Pathogenic CNVs', 'Non-pathogenic CNVs')) %>%
        ggplot(aes(id, n_overlap)) +
        geom_col(aes(fill = category), color = 'black', show.legend = FALSE) + 
        #geom_point() +
        # geom_smooth(aes(color = category, group = category))
        #geom_line(aes(color = category, group = category)) + 
        theme_minimal() + 
        facet_wrap(~category, nrow = 2)

    })
    

    
    
    output$df_overlap_cnvs <- renderDT({
      
      tmp_df <- df_overlap_cnvs_running() %>% filter(source == 'decipher')
      datatable(tmp_df, colnames = c('ID', 'Chrom', 'Start', 'End', 'Database', 'CNV size', 'Percentage Overlap (%)'),
                selection = 'single',
                options = list(
                  columnDefs = list(list(className = 'dt-center',  targets = 0:6))),  rownames= FALSE)
    })
    
    
    output$df_overlap_cnvs_nonpatho <- renderDT({
      
      tmp_df <- df_overlap_cnvs_running() %>% filter(source != 'decipher')
      datatable(tmp_df, colnames = c('ID', 'Chrom', 'Start', 'End', 'Database', 'CNV size', 'Percentage Overlap (%)'),
                selection = 'single',
                options = list(
                  columnDefs = list(list(className = 'dt-center',  targets = 0:6))),  rownames= FALSE)
    })
    
    upload_raw <- reactive({
      
      req(input$upload_bed_file)
      
      x <- input$upload_bed_file
      test231 <<-  input$filter_length
      
      tmp_path <-  x$datapath
      tmp_df <- read_tsv(tmp_path, col_names = c('chrom', 'start', 'end'))
      
      
    })
    
    
    running_upload <- reactive({
      
      req(upload_raw())

      test131311111 <<- upload_raw()
      tmp_df <- upload_raw() %>% mutate(size = end - start + 1) %>%
        filter(size >= input$filter_length)
      
      tmp_df

    })
    
    
    output$n_upload_cnv <- renderUI({
      
      req(running_upload())
      
      data_tmp <- upload_raw() %>% nrow()
      
      tablerStatCard(
        value =  data_tmp,
        title = "Number of uploaded CNVs",
        # trend = -10,
        width = 12
      )
      
      
    })
    
    output$n_upload_filter_cnv <- renderUI({
      
      req(upload_raw())
      req(running_upload())
      
      n_raw <- upload_raw() %>% nrow()
      n_filtered <- running_upload() %>% nrow()
      
      tablerStatCard(
        value =  paste0(n_filtered, '/', n_raw),
        title = "Number of filtered CNVs",
        # trend = -10,
        width = 12
      )
      
      
    })
    
    
    
    output$df_upload <- renderDataTable({
      
      req(running_upload())

      tmp_df <- running_upload() %>%
        mutate(score = pmap(list(chrom, start, end), function(a,b,c) get_model_score(a,b,c, hgcn_genes, model1))) %>%
        rowwise() %>%
        mutate(n_genes = score[[2]],
               score = score[[1]])
      
      test313313221 <<- tmp_df
      
      
      datatable(tmp_df, colnames = c('Chromosome', 'Start', 'End', 'Size', 'Score', 'nº genes'), rownames = FALSE)
      
      
    })
    
    
    
    running_de_novo <- reactive({
      
      req(input$start_analysis > 0)
      
      start_coordinates <- coord_user()[1]
      end_coordinates <- coord_user()[2]
      chrom_coordinates <- coord_user()[3]
      

      
      tmp_df <- denovo %>%
        filter(chrom == chrom_coordinates) %>%
        rowwise() %>%
        mutate(keep = Position %overlaps% c(start_coordinates, end_coordinates)) %>%
        ungroup() %>%
        filter(keep == TRUE) %>%
        select(-keep)
      
      tmp_df

    })
    
    
    
    output$df_de_novo <- renderDataTable({
      
    
      
      datatable(running_de_novo(), colnames = c('Chromosome', 'Position', 'Phenotype', 'Study name', 'PubmedID', 'Function Class'), rownames = FALSE)
      
      
    })
    
    
    output$ui_intersect <- renderUI({


      req(input$df_intersection_rows_selected)

      actionBttn(
        inputId = "take_intersect",
        label = "Run this region!",
        color = "warning",
        style = "material-flat",
        size = 'sm',
        block = TRUE
      )




    })

    
    output$download_dgenes <- downloadHandler(
      filename = function() {
        paste0('gene_panel', '_', Sys.Date(), ".tab")
      },
      content = function(file) {
        write.table(data_selected(), file, row.names = FALSE, sep = '\t')
      }
    )
    
    output$download_pubmed <- downloadHandler(
      filename = function() {
        paste0('pubmed_articles', '_', Sys.Date(), ".tab")
      },
      content = function(file) {
        write.table(running_pubmed(), file, row.names = FALSE, sep = '\t')
      }
    )
    
    output$download_de_novos <- downloadHandler(
      filename = function() {
        paste0('denovo_variants', '_', Sys.Date(), ".tab")
      },
      content = function(file) {
        write.table(running_de_novo(), file, row.names = FALSE, sep = '\t')
      }
    )
    
    output$download_df_overlap_cnvs <- downloadHandler(
      filename = function() {
        paste0('overlap_pathogenic_cnvs', '_', Sys.Date(), ".tab")
      },
      content = function(file) {
        write.table(df_overlap_cnvs_running_download_patho(), file, row.names = FALSE, sep = '\t')
      }
    )
    
    output$download_df_overlap_cnvs_nonpatho <- downloadHandler(
      filename = function() {
        paste0('overlap_non_pathogenic_cnvs', '_', Sys.Date(), ".tab")
      },
      content = function(file) {
        write.table(df_overlap_cnvs_running_download_nonpatho(), file, row.names = FALSE, sep = '\t')
      }
    )
    
    
    
  }
)