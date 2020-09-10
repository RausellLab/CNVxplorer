
tablerDashPage(
  
  navbar = tablerDashNav(
    
    shinyjs::useShinyjs(),
    useShinyalert(),
    
    
    
    navMenu = tablerNavMenu(
      
      tablerNavMenuItem(
        tabName = "overview",
        icon = "box",
        "Overview"
      ),
      tablerNavMenuItem(
        tabName = "genetic_evidence",
        icon = "box",
        "Genetic evidence"
      ),
      tablerNavMenuItem(
        tabName = "reg_region",
        icon = "box",
        "Regulatory elements"
      ),
      tablerNavMenuItem(
        tabName = "genomic_interactions",
        icon = "box",
        "Protein interaction network"
      ),
      tablerNavMenuItem(
        tabName = "disease",
        icon = "box",
        "Phenotypic analysis"
        
      ),
      tablerNavMenuItem(
        tabName = "model",
        icon = "box",
        "Model organism"
      ),


      tablerNavMenuItem(
        tabName = "pubmed",
        icon = "book",
        "Biomedical literature"
      ),
      tablerNavMenuItem(
        tabName = "tissue",
        icon = "box",
        "Tissue-specificity"
      ),
      
      tablerNavMenuItem(
        tabName = "fa",
        icon = "box",
        "Functional analysis"
      ),
      
      tablerNavMenuItem(
        tabName = "drugs",
        icon = 'book',
        "Druggable genome"
      ),
      tablerNavMenuItem(
        tabName = "docu",
        icon = "book",
        'Documentation'
      )
      
    ),
    id = "mymenu",
    src = "https://www.onlinelogomaker.com/applet_userdata/version2/5/0/18611424/projects/18611424.png",
    uiOutput('ref_user_genes_cnv'),
    uiOutput('counter_header'),
    uiOutput('ref_user_genes')
    
    # uiOutput('ref_user_filter_genes')
    
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
    copyrights = "@David Granjon, 2019"),
  
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
        shinyjs::useShinyjs(),
        use_waiter(),
        setShadow(class = 'card'),
        chooseSliderSkin("Nice"),
        
        fluidRow(
          column(
            width = 3,

            tablerCard(width = 12,
                       title = NULL,
                       collapsible = FALSE,
                       closable = FALSE,
                       zoomable = FALSE,
                       statusSide = 'left',
                       prettyRadioButtons(
                         inputId = "input_geno_karyo",
                         label = "Choose:", 
                         choices = c("Genomic coordinates", "G banding", 'Multiple coordinates'),
                         inline = TRUE, 
                         status = "primary",
                         fill = TRUE
                       ),
                       conditionalPanel(
                         condition = "input.input_geno_karyo != 'Multiple coordinates'",
                       selectizeInput(inputId = 'input_chrom', label = 'Chromosome', choices = human_chrom,
                                      selected = NULL, multiple = FALSE,
                                      options = NULL),
                       uiOutput('choose_geno_karyo1'),
                       uiOutput('choose_geno_karyo2')),
                       conditionalPanel(
                         condition = "input.input_geno_karyo == 'Multiple coordinates'",
                         fileInput("file_cnv", label = h5("Upload file:")),
                         downloadLink('download_file_1', label = "Download file example [#1]")
                       ),
                       
                       
            ),

            actionBttn(
              inputId = "start_analysis",
              label = "Run!",
              color = "success",
              style = "material-flat",
              size = 'lg',
              block = TRUE
            ),

            tags$br(),
            tags$br(),
            tags$br()
          ),
          column(6,
                 conditionalPanel(
                   condition = "input.input_geno_karyo == 'Multiple coordinates'",
                 tablerCard(width = 12,
                            title = 'Input file',
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            statusSide = 'left',


                     DTOutput('cnv_file')),
   
                 ),
                 plotOutput('plot_chrom', height = 200)

                 
                 
                 
                 
          ),
          
          column(
            width = 3,
            
            conditionalPanel(
              condition = "input.input_geno_karyo == 'Multiple coordinates'",

              tablerCard(width = 12,
                         title = NULL,
                         collapsible = FALSE,
                         closable = FALSE,
                         zoomable = FALSE,
                         statusSide = 'left',
                         selectizeInput(inputId = 'select_all_cnvs',
                                        label = tags$b('Select all the variants?'),
                                        choices = c('Yes' = 'yes', 'No' = 'no'),
                                        selected = NULL,
                                        multiple = FALSE,
                                        options = NULL),
           uiOutput('n_variants'),
           uiOutput('ref_user_region_file'),
           uiOutput('checking_quality_file')
                         
              ), 
              ),
            tablerCard(width = 12,
                       title = NULL,
                       collapsible = FALSE,
                       closable = FALSE,
                       zoomable = FALSE,
                       statusSide = 'left',
                       pickerInput(
                         inputId = "11type_cnv",
                         label = tags$b("CNV type"), 
                         choices = list("-" = NA, "Deletion" = 1, "Duplication" = 0)
                       ),
                       pickerInput(
                         inputId = "11denovo_yes_no",
                         label = tags$b("de novo CNV?"), 
                         choices = list("-" = NA, 'Yes' = 1, 'No' = 0)
                       ),
                       
            ), 
           
            uiOutput('ref_user_region'),
            uiOutput('ref_user_cytoband'),
            uiOutput('ref_user_length'),
            uiOutput('check_blacklist')

          )
        ),
        
        tags$hr(),
        fluidRow(
          column(width = 6,
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'success',
                            statusSide = 'left',
                            title = tagList(shiny::icon("database"), "Overlap with reference CNV databases"),
                            fluidRow(
                              column(width = 4,
                                     uiOutput('n_syndromes')),
                              
                              column(width = 4,
                                     uiOutput('n_cnv_patho')),
                              column(width = 4,
                                     uiOutput('n_cnv_nopatho'))))
          ),
          
          column(width = 6,
                 
                 
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'primary',
                            statusSide = 'left',
                            title = tagList(shiny::icon("database"), "Overlap with disease associated genes"),
                            fluidRow(
                              column(width = 4,
                                     uiOutput('n_disease')),
                              
                              column(width = 4,
                                     uiOutput('n_clinvar')),
                              column(width = 4,
                                     uiOutput('n_gwas')
                              ))))),
        
        
        fluidRow(
          
          column(width = 6,
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'success',
                            statusSide = 'left',
                            title = tagList(shiny::icon("book"), "Mouse model information"),
                            fluidRow(
                              column(width = 6,
                                     
                                     uiOutput('n_mortality')),
                              column(width = 6,
                                     uiOutput('n_embryo')
                              ))),
                 fluidRow(
                   tablerCard(width = 12,
                              collapsible = FALSE,
                              closable = FALSE,
                              zoomable = FALSE,
                              status = 'success',
                              statusSide = 'left',
                              title = tagList(shiny::icon("hospital"), "Clinical information"), 
                              fluidRow(
                                column(width = 6,
                                       
                                       uiOutput('hpo_unique')), 
                                column(width = 6,
                                       uiOutput('hpo_unique_genes') 
                                ))))

                 
          ),
          column(width = 6,
                 
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'success',
                            statusSide = 'left',
                            title = tagList(shiny::icon("book"), "Regulatory elements disrupted"),
                            
                            
                            
                            fluidRow(
                              
                              column(width = 4,
                                     uiOutput('n_enhancer'),
                                     uiOutput('n_tfs')),
                              
                              column(width = 4,
                                     uiOutput('n_tads')),
                              
                              column(width = 4,
                                     uiOutput('n_lncrna'),
                                     uiOutput('n_mirna')
                              )
                              
                              
                            )),
                 
                 
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'success',
                            statusSide = 'left',
                            title = tagList(shiny::icon("book"), "Biomedical literature"),
                            
                            
                            
                            fluidRow(
                              
                              column(width = 6,
                                     uiOutput('n_pubmed_del')),
                              column(width = 6,
                                     uiOutput('n_pubmed_dup')
                              )
                              
                              
                            ))
                 
          )
        ),
        
        
        fluidRow(
          
          
          
          column(width = 12,
                 
                 tags$hr(),
                 
                 tablerCard(
                   title = "Comparison CNV size with other CNVs databases",
                   plotOutput('plot_size'),
                   width = 12,
                   overflow = TRUE,
                   collapsible = FALSE,
                   closable = FALSE,
                   options = tagList(
                     
                     prettyRadioButtons(
                       inputId = "select_density",
                       label = '', 
                       choices =  split(c("local", "global"), c('Local', 'Global')),
                       selected = 'local',
                       inline = TRUE, 
                       status = "primary",
                       fill = TRUE
                     )
                     
                     
                   )
                 )
          )
        ),
        fluidRow(
          
          
          column(width = 12,
                 
                 fluidRow(
                   column(width = 6,
                          uiOutput('n_enhdancer'),
                          uiOutput('n_lnscrna')),
                   column(width = 6,
                          uiOutput('n_tdads')
                   ))
          )
          
        )
        
      ),
      
      tablerTabItem(
        tabName = "fa",

        tablerCard(
          title = "Functional Profile",
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          uiOutput('plot_df') ,
          options = tagList(
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
        tabName = "genetic_evidence",
        tablerCard(
          title = "Overlap with CNV Syndromes",
          DTOutput("cnv_syndromes"),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            uiOutput('ui_select_cnv_syndrome')
            # downloadButton("downl323rtfgsoad_dgenes", "Download table")
          )
        ),
        tablerCard(
          title = "Overlap with pathogenic/likely pathogenic CNVs (DECIPHER)",
          DTOutput('df_overlap_cnvs'),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            # downloadButton("download_df_overlap_cnvs", "Download table")
            
          )
        ),
        tablerCard(
          title = "Overlap with non-pathogenic CNVs",
          DTOutput('df_overlap_cnvs_nonpatho'),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            uiOutput('select_db_no_patho_cnvs')
            # downloadButton("download_df_overlap_cnvs_nonpatho", "Download table")
            
          )
        ),
        tablerCard(
          title = "Intersection of disease gene databases",
          plotOutput("plot_upset_disease"),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE
        ),

        fluidRow(
          tablerCard(
            title = "Overlap with disease genes",
            DTOutput("dgenes") %>% withSpinner(type = 5),
            width = 4,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE
          ),
          tablerCard(
            title = "Disease evidences",
            DTOutput("select_gene_disease") %>% withSpinner(type = 5),
            width = 8,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('input_source')
              # downloadButton("downloadfdsfs_dgenes", "Download table")
            )
          )),
        tablerCard(
          title = "Overlap with non-disease genes",
          DTOutput("dgenes_no_disease") %>% withSpinner(type = 5),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            # downloadButton("download_adsasdasdasdasdasdgenes", "Download table")
          )
        ),


        fluidRow(

          tablerCard(
            title = "Disease-associated variants",
            DTOutput('df_variants'),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('ui_select_clinvar_gwas')
              # downloadButton("downlda323rtfgsoad_dgenes", "Download table")
            )
          ),
          tablerCard(
            title = "De novo variants (denovo-db)",
            DTOutput('df_de_novo'),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              # downloadButton("download_dadasdasdasdasdasdasde_novos", "Download table")
              
            )
          )

          
        )
      ),
      tablerTabItem(
        tabName = "reg_region",
        fluidRow(

          
          tablerCard(
            title = "Intersection of disease target-genes from regulatory elements",
            plotOutput("plot_upset_disease_reg"),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE
          ),
          column(12,
          fluidRow(
            tablerCard(
              title = "Overlap with disease target-genes from regulatory elements",
              DTOutput("dgenes_reg") %>% withSpinner(type = 5),
              width = 4,
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
                # downloadButton("download_dgenes", "Download table")
              )
            ),
            tablerCard(
              title = "Disease evidences",
              DTOutput("select_gene_disease_reg") %>% withSpinner(type = 5),
              width = 8,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE,
              options = tagList(
                uiOutput('input_source_reg')
                # downloadButton("downloadfdsfs_dgenes", "Download table")
              )
            ))
          ),
          
          tablerCard(
            title = "Non-disease target genes from disrupted regulatory elements",
            DTOutput("genes_from_reg_regions") %>% withSpinner(type = 5),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('choose_reg_region')
            )
          ),
          tablerCard(title = 'Enhancers disrupted',
                     collapsible = FALSE,
                     closable = FALSE,
                     dataTableOutput('df_enhancer') %>% withSpinner(type = 5),
                     width = 12,
                     options = tagList(
                       uiOutput('n_target_enh_notmap'),
                       uiOutput('n_target_enh'),
                       uiOutput('switch_enhancers')
                       
                     ))),
        fluidRow(width = 12, 
                 tablerCard(title = 'Conservation - Phast100way histogram',
                            plotOutput('p100_enhancer'),
                            collapsible = FALSE,
                            closable = FALSE,
                            width = 4),
                 tablerCard(title = 'Conservation - Phast46way placental histogram',
                            plotOutput('p46pla_enhancer'),
                            collapsible = FALSE,
                            closable = FALSE,
                            width = 4),
                 tablerCard(title = 'Conservation - Phast46way primates histogram',
                            plotOutput('p46pri_enhancer'),
                            collapsible = FALSE,
                            closable = FALSE,
                            width = 4)
        ),
        

        tablerCard(title = 'micro-RNAs (miRNAs) disrupted',
                   DTOutput('df_mirna'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_mirna_notmap'),
                     uiOutput('n_target_mirna'),
                     uiOutput('switch_mirnas')
                   )),
        tablerCard(title = 'Transcription factors (TFs) disrupted',
                   DTOutput('tf_df'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_tf_notmap'),
                     uiOutput('n_target_tf'),
                     uiOutput('switch_tfs')
                   )),

        tablerCard(title = 'Long noncoding RNAs (lncRNAs) disrupted',
                   DTOutput('lncrna_df'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_lncrna_notmap'),
                     uiOutput('n_target_lncrna'),
                     uiOutput('switch_lncrnas')
                   )),
        tablerCard(title = 'Topologically Associating Domains (TADs) disrupted',
                   DTOutput('df_tads'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('ui_tad')
                     # uiOutput('switch_tads')
                   )),
        
      ),
      
      tablerTabItem(
        tabName = "genomic_interactions",
        
        # tagAppendAttributes(actionButton(ns("open"), "Network Options"), class = "btn-outline-primary"),
        

        tablerCard(title = 'Network',
                   width = 12,
                   forceNetworkOutput("network_ppi")
        )
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
      
      # tablerTabItem(
      #   tabName = "gene",
      #   fluidRow(
      #     tablerCard(title = '19,281 protein-coding genes ',
      #                reactableOutput("table_all_genes"),
      #                width = 12)),
      #   tablerCard(title = 'RNA Expression (GTEx)',
      #              DTOutput('lncrn2a_df'),
      #              width = 12)
      #   
      # ),
      
      tablerTabItem(
        tabName = "tissue",
        fluidRow(
          column(width = 3,
                 tablerCard(title = 'Filter options:',
                            tags$b('Display by:'),
                            selectInput(
                              "gtex_gene_tissue", "",
                              c(Gene = "Gene",
                                Tissue = "Tissue"),
                              selected = 'Tissue'),
                            conditionalPanel(
                              condition = "input.gtex_gene_tissue == 'Gene'",
                              uiOutput('gtex_gene')),
                            conditionalPanel(
                              condition = "input.gtex_gene_tissue == 'Tissue'",
                              uiOutput('gtex_tissue')),
                            
                            width = 12)
                 
                 
          ),
          
          tablerCard(title = 'RNA Expression (GTEx)',
                     plotOutput('tissue_gtex'),
                     width = 9)
          
          
          
        ),

        
        
        

        fluidRow(
          column(width = 3,
                 tablerCard(title = 'Configuration:',
                            tags$b('Reference database:'),
                            selectInput(
                              "tissue_expression_dbs", "",
                              c("GTEx" = "gtex",
                                "ENCODE mouse" = "encode_mouse")),
                            width = 12)
                 
                 
          ),
          
          tablerCard(title = 'Tissue-specific enrichment analysis (TSEA)',
                     plotOutput('plot_tsea'),
                     width = 9,
                     options = tagList(
                       
                       
                       switchInput(
                         inputId = "enable_tsea",
                         label = "Run?",
                         value = FALSE,
                         onStatus = "success",
                         offStatus = "danger",
                         width = 'auto',
                         size = 'mini'
                         
                       )
                       
                     ))
          
          
        ),
        fluidRow(
          column(width = 3,
                 tablerCard(title = 'Filter options:',
                            tags$b('Filter by gene:'),
                            selectInput(
                              "gene_yes_no", "",
                              c(No = "No",
                                Yes = "Yes")),
                            conditionalPanel(
                              condition = "input.gene_yes_no == 'Yes'",
                              uiOutput('gene_tissue')),
                            tags$b('Filter by tissue:'),
                            selectInput(
                              "tissue_yes_no", "",
                              c(No = "No",
                                Yes = "Yes")),
                            conditionalPanel(
                              condition = "input.tissue_yes_no == 'Yes'",
                              uiOutput('select_tissue')),
                            
                            width = 12)
                 
          ),
          tablerCard(title = 'Protein Expression (Human Protein Atlas)',
                     DTOutput('tissue_hpa'),
                     width = 9))
        
      ),
      tablerTabItem(
        tabName = "disease",
        fluidRow(
          uiOutput('n_hp_chosen')

          
        ),
        fluidRow(
          tablerCard(title = 'Phenotype terms',
                     width = 8,
                     closable = FALSE,
                     collapsible = FALSE,
                     zoomable = FALSE,
                     multiInput(
                       inputId = "chosen_hp",
                       label = "",
                       choices = NULL,
                       width = '100%',
                       choiceNames = vector_total_terms$term_desc,
                       choiceValues = vector_total_terms$term
                     ),
                     options = tagList(
                       actionBttn(
                         inputId = "reset_pheno_analysis",
                         label = "Reset",
                         size = 'sm',
                         color = "default",
                         style = "material-flat",
                         block = TRUE
                       )
                     )
                     
          ), 
          column(width = 4,
                 
                 tablerCard(title = 'Find HP terms in a clinical text',
                            width = 12,
                            closable = FALSE,
                            collapsible = FALSE,
                            zoomable = FALSE,
                           textInput(inputId = 'text_recognition', label = NULL, value = "", width = NULL,
                                     placeholder = 'Enter your clinical annotations...'),
                           tags$hr(),
                           DTOutput('entities_df')

                           ),
                 

                          
                            

                 
                 tablerCard(title = 'Suggested phenotype terms',
                            width = 12,
                            DTOutput('suggest_df')
                            
                 )),

          column(width = 12,
                 fluidRow(width = 12,

                          tablerCard(title = 'Gene-disease associations',
                                     closable = FALSE,
                                     collapsible = FALSE,
                                     zoomable = FALSE,
                                     DTOutput('dt_running_sim_score') %>% withSpinner(type = 5),
                                     width = 12,
                                     options = tagList(
                                       uiOutput('n_genes'),
                                       uiOutput('n_diseases'),
                                       
                                       selectizeInput(inputId = 'input_inheritance',
                                                      label = tags$b('Filter by mode of inheritance:'),
                                                      choices = vector_inheritance,
                                                      selected = NULL,
                                                      multiple = FALSE,
                                                      options = NULL)
                                       
                                       
                                     )),
                          # fluidRow(width = 12,
                                   
                          tablerCard(title = 'HPO terms associated with genes',
                                     DTOutput('hpo_assoc_genes'),
                                     width = 6),
                          tablerCard(title = 'HPO terms associated with diseases',
                                     DTOutput('hpo_assoc_diseases'),
                                     width = 6)
                          # )
                 )),
          
          tablerCard(title = 'Anatomical entities associated with HPO terms',
                     highchartOutput('plot_anatomy'),
                     width = 12),
          tablerCard(title = 'Phenotypic similarity score',
                     
                     plotOutput('plot_similarity_genes'),
                     width = 12,
                     options = tagList( 
                       
                       prettyRadioButtons(
                         inputId = "select_sim_gene_disease",
                         label = '', 
                         choices = list('Genes' = 'genes', 'OMIM diseases' = 'diseases'),
                         inline = TRUE, 
                         status = "primary",
                         fill = TRUE
                       ))
          )

        )
      ),
      tablerTabItem(
        tabName = "model",
        tablerCard(title = 'Mouse phenotypes associated with genes',
                   plotOutput('agg_model'),
                   width = 12),
        fluidRow(

          tablerCard(title = 'Genes associated with phenotypes in mouse (MGI)',
                     DTOutput('model_genes'),
                     width = 12))
        
      ),
      tablerTabItem(
        tabName = "drugs",
        fluidRow(
          tablerCard(title = 'Approved drugs associated with genes',
                     DTOutput('drugbank_df'),
                     width = 12)

          )
        
      ),
      tablerTabItem(
        tabName = "docu",
      
          tablerCard(
            collapsible = FALSE,
            closable = FALSE,
            title = 'Overview',
            # tags$img(src='overview.jpg'),
            includeMarkdown('doc/documentation.Rmd'),
                     width = 12)
        
      ),
      tablerTabItem(
        tabName = "pubmed",
        fluidRow(
          tablerCard(
            title = "Pubmed articles associated with the region",
            
            collapsible = FALSE,
            closable = FALSE,
            DTOutput('del_dup_pubmed') %>% withSpinner(type = 5),
            overflow = TRUE,
            width = 12,
            options = tagList(
              uiOutput('ui_only_omim'),
              prettyRadioButtons(
                inputId = "select_del_dup",
                label = '',   
                choices =  split(c('deletions', 'duplications'), c('Deletions', 'Duplications')),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              )
            )
          ),
          tablerCard(
            title = "Abstract",
            
            htmlOutput("abstract_html"),
            width = 8,
            overflow = TRUE
          ),
          tablerCard(
            title = "Entities found",
            
            DTOutput("abstract_df"),
            width = 4,
            overflow = TRUE
          ),
          tablerCard(
            collapsible = FALSE,
            closable = FALSE,
            title = NULL,
            width = 3,
            sliderInput("min_threshold_cooccurrence", label = tags$b("Select minimum co-occurrence:"), min = 2, 
                        max = 50, value = 2)
          ),
            
          tablerCard(
            collapsible = FALSE,
            closable = FALSE,
            title = "Co-occurrence network",
            
            
            plotOutput("plot_net_pubmed") %>% withSpinner(type = 5),
            width = 9,
            overflow = TRUE,
            options = tagList(

              prettyRadioButtons(
                inputId = "choose_net_source",
                label = '',
                choices = list('Both' = 'both', 'Deletions' = 'deletion', 'Duplications' = 'duplication'),
                # selected = 'total',
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              ),
              switchInput(
                inputId = "enable_net",
                label = "Run?",
                value = FALSE,
                onStatus = "success",
                offStatus = "danger",
                width = 'auto',
                size = 'mini'
                
              )
          ))

          
        )
        
      ),
      tablerTabItem(
        tabName = "down_report",
        fluidRow(
          tablerCard(
            title = "Personalize your report:",
            
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
)