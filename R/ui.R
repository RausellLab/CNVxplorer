
tablerDashPage(
  
  navbar = tablerDashNav(
    img(src = 'imagine_logo.png', height = '75px', width = '150px'),
    shinyjs::useShinyjs(),
    useShinyalert(),
    tags$head(includeHTML("www/google_analytics.html")),

    navMenu = tablerNavMenu(
      
      tablerNavMenuItem(
        tabName = "overview",
        icon = "box",
        "Overview"
      ),
      tablerNavMenuItem(
        tabName = "genetic_evidence",
        icon = "box",
        "Clinical genetics evidence"
      ),
      tablerNavMenuItem(
        tabName = "reg_region",
        icon = "box",
        "Regulatory regions"
      ),
      tablerNavMenuItem(
        tabName = "disease",
        icon = "box",
        "Phenotypic analysis"
        
      ),

      tablerNavMenuItem(
        tabName = "model",
        icon = "box",
        "KO mouse phenotypes"
      ),

      tablerNavMenuItem(
        tabName = "fa",
        icon = "box",
        "Functional enrichment analysis"
      ),
      tablerNavMenuItem(
        tabName = "tissue",
        icon = "box",
        "Tissue expression patterns"
      ),
      
      tablerNavMenuItem(
        tabName = "genomic_interactions",
        icon = "box",
        "Protein interaction"
      ),
      tablerNavMenuItem(
        tabName = "pubmed",
        icon = "book",
        "Biomedical literature"
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
  ),
  footer = tablerDashFooter(
    tablerIcon(name = "maestro", lib = "payment"),
    tablerIcon(name = "mastercard", lib = "payment"),
    copyrights = "@David Granjon, 2019"),
  
  title = "CNVxplorer - A web tool for the clinical interpretation of CNVs",
  body = tablerDashBody(
    
    # link js
    tags$head(tags$link(includeScript("www/index.js"))),
    tags$head(tags$style("a{cursor:pointer;}")),
    
    
    tags$head(tags$script('
  $(document).on("shiny:sessioninitialized", function(event) {
    $(\'a[data-value="Page1"]\').tab("show");
  });
')),
    
    tablerTabItems(
      
      
      
      tablerTabItem(
        tabName = "overview",
        shinyjs::useShinyjs(),
        # use_waiter(),
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
                         label = "Choose: (GRCh37 coordinates)", 
                         choices = c("Genomic coordinates", "G banding", 'Multiple coordinates (NGS)'),
                         inline = TRUE, 
                         status = "primary",
                         fill = TRUE
                       ),
                       conditionalPanel(
                         condition = "input.input_geno_karyo != 'Multiple coordinates (NGS)'",
                       selectizeInput(inputId = 'input_chrom', label = 'Chromosome', choices = human_chrom,
                                      selected = NULL, multiple = FALSE,
                                      options = NULL),
                       uiOutput('choose_geno_karyo1'),
                       uiOutput('choose_geno_karyo2')),
                       conditionalPanel(
                         condition = "input.input_geno_karyo == 'Multiple coordinates (NGS)'",
                         fileInput("file_cnv", label = h5("Upload file (.tsv):")%>% helper(type = "inline",
                                                                                                               # icon = "exclamation",
                                                                                                               style = "text-indent: 0.5em;",
                                                                                                               title = "File missing",
                                                                                                               size = "m",
                                                                                                               buttonLabel = 'OK',
                                                                                                               content = c("If you can not find your file, please make sure the file extension is .bed or .tsv")
                         ), accept = c('.bed', '.tsv')),
                         downloadLink('download_file_1', label = "Download file example [#1]")
                       )
                       
                       
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
                   condition = "input.input_geno_karyo == 'Multiple coordinates (NGS)'",
                 tablerCard(width = 12,
                            title = 'Input variants',
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            statusSide = 'left',


                     DTOutput('cnv_file'))
   
                 ),
                 plotOutput('plot_chrom', height = 200),
                 DTOutput('cnvscore_rules')

          ),
          
          column(
            width = 3,
            
            conditionalPanel(
              condition = "input.input_geno_karyo == 'Multiple coordinates (NGS)'",

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
                         conditionalPanel(
                           condition = "input.select_all_cnvs == 'no'",
                           tags$p(tags$b('To refresh the app after each row selection, you need to click on the button "Run".'))),
           uiOutput('n_variants'),
           uiOutput('ref_user_region_file'),
           uiOutput('checking_quality_file')
           
                         
              )
              ),
            uiOutput('ref_user_cnvscore'),
            # uiOutput('cnvscore_rules'),
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
                            title = tagList(shiny::icon("database"),
                                            HTML("<a onclick=","customHref('genetic_evidence')", 'target="_top"' ,">", 
                                                 "Overlap with reference CNV databases (Link)","</a>")),
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
                            title = tagList(shiny::icon("database"),
                                            HTML("<a onclick=","customHref('genetic_evidence')", 'target="_top"' ,">", 
                                                 "Overlap with disease associated genes and SNV variants (Link)","</a>")),
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
                            title = tagList(shiny::icon("book"),
                                            HTML("<a onclick=","customHref('model')", 'target="_top"' ,">", 
                                                 "Mouse model information (Link)","</a>")),
                            fluidRow(
                              column(width = 6,
                                     
                                     uiOutput('n_mortality')),
                              column(width = 6,
                                     uiOutput('n_embryo')
                              ))),
                   tablerCard(width = 12,
                              collapsible = FALSE,
                              closable = FALSE,
                              zoomable = FALSE,
                              status = 'success',
                              statusSide = 'left',
                              title = tagList(shiny::icon("hospital"),
                                              HTML("<a onclick=","customHref('disease')", 'target="_top"' ,">", 
                                                   "Clinical information (Link)","</a>")),
                              fluidRow(
                                column(width = 6,
                                       
                                       uiOutput('hpo_unique')), 
                                column(width = 6,
                                       uiOutput('hpo_unique_genes') 
                                )))

                 
          ),
          column(width = 6,
                 
                 tablerCard(width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            zoomable = FALSE,
                            status = 'success',
                            statusSide = 'left',
                            title = tagList(shiny::icon("book"),
                                    HTML("<a onclick=","customHref('reg_region')", 'target="_top"' ,">", 
                                         "Overlap with regulatory elements and TADs (Link)","</a>")),
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
                            title = tagList(shiny::icon("book"),
                                            HTML("<a onclick=","customHref('pubmed')", 'target="_top"' ,">", 
                                                 "Biomedical literature (Link)","</a>")),
                            
                            
                            
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
                   title = "Comparison CNV size with other CNVs databases" %>% helper(type = "inline",
                                                                                      style = "text-indent: 0.5em;",
                                                                                      title = "Two options: global & local",
                                                                                      size = "m",
                                                                                      buttonLabel = 'OK',
                                                                                      content = c("
CNVxplorer compares the length of the CNV provided by the user and the length distribution of the CNVs found in five databases (DGV, gnomAD, DECIPHER Control, DECIPHER, ClinVar). To make a comparison, we provide two approaches:",
                                                                                                  
                                                                                                  "<b> Global: </b> The length of the CNV is compared with the total number of CNVs available in the databases.",
                                                                                                  "<b> Local: </b> The comparison is made exclusively with the CNVs mapping the query.",
                                                                                                  "Both options usually provide similar results. But this does not always have to be the case, especially on those queries with a short length.")
                                                                                      
                   ),
                   plotOutput('plot_size') ,
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
          title = "GO term annotation" %>% helper(type = "inline",
                                                  style = "text-indent: 0.5em;",
                                                  title = "Functional profile",
                                                  size = "m",
                                                  buttonLabel = 'OK',
                                                  content = c("This panel does not provide a functional enrichment analysis but the functional annotation of the selected genes. The specificity of the annotation can be increased by setting the parameter 'Level'")
                                                  ),
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
          title = "Overlap with CNV Syndromes" %>% helper(type = "inline",
                                                          style = "text-indent: 0.5em;",
                                                          title = "CNV Syndromes",
                                                          size = "m",
                                                          buttonLabel = 'OK',
                                                          content = 
                                                            c("<b>How is the overlap calculated?</b>",
                                                               "To calculate the overlap, we use CNV Syndromes from ClinGen and DECIPHER  as a reference and, as a query, the CNV(s) entered by the user. For instance, a 100% overlap means that the CNV syndrome is completely mapping the user's CNV. ",
                                                               "In the Documentation - FAQ tab, you can find a picture where we illustrate this point",
                                                               "<b>Duplicated entries</b>",
                                                               "We collect CNV Syndromes from two curated databases: ClinGen and DECIPHER. These two databases share multiple entries. Therefore, it is expected that you will find some entries in both sources.")

                                                            ),
          DTOutput("cnv_syndromes"),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            uiOutput('ui_select_cnv_syndrome')
          )
        ),
        tablerCard(
          title = "Overlap with pathogenic/likely pathogenic CNVs (DECIPHER & ClinVar)" %>% helper(type = "inline",
                                                                                     style = "text-indent: 0.5em;",
                                                                                     title = "Likely pathogenic/pathogenic CNVs",
                                                                                     size = "m",
                                                                                     buttonLabel = 'OK',
                                                                                     content = c("<b> Overlap calculation </b>",
                                                                                                  "To calculate the overlap, we use likely pathogenic/pathogenic CNVs from DECIPHER and ClinVar as a reference and, as a query, the CNV(s) entered by the user. For instance, a 100% overlap means that the DECIPHER's CNV is completely mapping the user's CNV. ",
                                                                                                  "In the Documentation - FAQ tab, you can find a picture where we illustrate this point.",
                                                                                                  "<b> Phenotypic similarity </b>",
                                                                                                  "For every CNV associated with HP terms, you can calculate the phenotypic similarity with the patient’s clinical symptoms. To do so, go to the Phenotypic analysis tab, enter HP terms and you will find a panel with the list of DECIPHER CNVs and their respective scores.")
                                                                                     
          ),
          DTOutput('df_overlap_cnvs'),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            uiOutput('ui_select_decipher_clinvar')
          )
        ),
        tablerCard(
          title = "Overlap with non-pathogenic CNVs"  %>% helper(type = "inline",
                                                                 style = "text-indent: 0.5em;",
                                                                 title = "Likely pathogenic/pathogenic CNVs",
                                                                 size = "m",
                                                                 buttonLabel = 'OK',
                                                                 content = c("To calculate the overlap, we use the user's CNV(s) as a reference and, as a query, the non-pathogenic CNVs from databases. For instance, a 100% overlap means that the user's query is completely mapping the non-pathogenic CNV.",
                                                                             "In Documentation - FAQs, you can find a picture where we illustrate this point.")
          ),
                                                                 
          DTOutput('df_overlap_cnvs_nonpatho'),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE,
          options = tagList(
            uiOutput('select_db_no_patho_cnvs')
          )
        ),
        tablerCard(
          title = "Overlap with disease-associated genes from different databases" %>% helper(type = "inline",
                                                                      style = "text-indent: 0.5em;",
                                                                      title = "Intersection of disease evidence",
                                                                      size = "m",
                                                                      buttonLabel = 'OK',
                                                                      content = c("To assess whether a gene is disease-associated, we collect information from 5 databases. This panel shows, for each gene, the number of databases that support this association. The maximum number /(five/) represents that the gene-disease association is well supported and widely known.")
          ),
          plotOutput("plot_upset_disease"),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE
        ),
          tablerCard(
            title = "Overlap with disease genes",
            DTOutput("dgenes"),
            width = 4,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE
          ),
          tablerCard(
            title = "Disease evidence",
            DTOutput("select_gene_disease") %>% withSpinner(type = 5),
            width = 8,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('input_source')
            )
          ),
        tablerCard(
          title = "Overlap with disease & non-disease genes"  %>% helper(type = "inline",
                                                                         style = "text-indent: 0.5em;",
                                                                         title = "Pathogenicity scores",
                                                                         size = "m",
                                                                         buttonLabel = 'OK',
                                                                         content = c("To facilitate the interpretation of the pathogenicity prediction scores, we transform each one of them into a percentile scale. For each score, percentiles near 100 reflect greater pathogenicity.",
                                                                                     "You can find more information about each score in Documentation - FAQs")),
                                                                         
          DTOutput("dgenes_no_disease") %>% withSpinner(type = 5),
          width = 12,
          collapsible = FALSE,
          closable = FALSE,
          overflow = TRUE
        ),


          tablerCard(
            title = "Disease-associated variants",
            DTOutput('df_variants'),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('ui_select_clinvar_gwas')
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

            )
          )

      ),
      tablerTabItem(
        tabName = "reg_region",
        fluidRow(

          
          tablerCard(
            title = "Intersection of disease target-genes databases (from overlapping regulatory elements)",
            plotOutput("plot_upset_disease_reg"),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE
          ),
          column(12,
          fluidRow(
            tablerCard(
              title = "Disease target-genes from regulatory elements",
              DTOutput("dgenes_reg") %>% withSpinner(type = 5),
              width = 4,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE
            ),
            tablerCard(
              title = "Disease evidence" %>% helper(type = "inline",
                                                     style = "text-indent: 0.5em;",
                                                     title = "Intersection of disease evidence",
                                                     size = "m",
                                                     buttonLabel = 'OK',
                                                     content = c("To assess whether a gene is disease-associated, we collect information from 5 databases. This panel shows, for each gene, the number of databases that support this association. The maximum number (five) represents that the gene-disease association is well supported and widely known.")),
              DTOutput("select_gene_disease_reg") %>% withSpinner(type = 5),
              width = 8,
              collapsible = FALSE,
              closable = FALSE,
              overflow = TRUE,
              options = tagList(
                uiOutput('input_source_reg')
              )
            ))
          ),
          
          tablerCard(
            title = "Disease & non-disease target genes from overlapping regulatory elements" %>% helper(type = "inline",
                                                                                                       style = "text-indent: 0.5em;",
                                                                                                       title = "Pathogenicity scores",
                                                                                                       size = "m",
                                                                                                       buttonLabel = 'OK',
                                                                                                       content = c("To facilitate the interpretation of the pathogenicity prediction scores, we transform each one of them into a percentile scale. For each score, percentiles near 100 reflect greater pathogenicity.",
                                                                                                                   "You can find more information about each score in Documentation - FAQs")),
            DTOutput("genes_from_reg_regions") %>% withSpinner(type = 5),
            width = 12,
            collapsible = FALSE,
            closable = FALSE,
            overflow = TRUE,
            options = tagList(
              uiOutput('choose_reg_region')
            )
          ),
          tablerCard(title = 'Overlapping enhancers',
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
        

        tablerCard(title = 'Overlapping micro-RNAs (miRNAs)',
                   DTOutput('df_mirna'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_mirna_notmap'),
                     uiOutput('n_target_mirna'),
                     uiOutput('switch_mirnas')
                   )),
        tablerCard(title = 'Overlapping Transcription factors (TFs)',
                   DTOutput('tf_df'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_tf_notmap'),
                     uiOutput('n_target_tf'),
                     uiOutput('switch_tfs')
                   )),

        tablerCard(title = 'Overlapping long noncoding RNAs (lncRNAs)',
                   DTOutput('lncrna_df'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('n_target_lncrna_notmap'),
                     uiOutput('n_target_lncrna'),
                     uiOutput('switch_lncrnas')
                   )),
        tablerCard(title = 'Overlapping Topologically Associating Domains (TADs)',
                   DTOutput('df_tads'),
                   collapsible = FALSE,
                   closable = FALSE,
                   width = 12,
                   options = tagList(
                     uiOutput('switch_tads')
                   ))
        
      ),
      
      tablerTabItem(
        tabName = "genomic_interactions",
        column(width = 9,
        tablerCard(title = 'Protein interaction network' %>% helper(type = "inline",
                                                                    style = "text-indent: 0.5em;",
                                                                    title = "Interactions evidence",
                                                                    size = "m",
                                                                    buttonLabel = 'OK',
                                                                    content = c("CNVxplorer represents only interactions from STRING with a high evidence score (equal or higher than 700)")),
                   width = 12,
                   forceNetworkOutput("network_ppi")
        )),
        column(width = 3,
               tablerCard(title = 'Filter by gene', width = 12,
                          prettyRadioButtons(
                            inputId = "filter_by_gene_ppi",
                            label = "Choose:", 
                            choices = c("No", "Yes"),
                            inline = TRUE, 
                            status = "primary",
                            fill = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.filter_by_gene_ppi == 'Yes'",
                            uiOutput("output_select_gene"))
               ),
               tablerCard(title = 'Color legend', width = 12,
                          htmlOutput('legend_html'))
               
),
               
        tablerCard(title = 'Nº of protein-protein interactions',
                   width = 12,
                   plotOutput('frequency_network')
        )
      ),
          
          
      
      
      tablerTabItem(
        tabName = "cnv_ngs",
        fluidRow(
          column(4,
                 tablerCard(title = 'Input variants' %>% helper(type = "inline",
                                                            # icon = "exclamation",
                                                            style = "text-indent: 0.5em;",
                                                            title = "+1 start - genome interval",
                                                            size = "m",
                                                            buttonLabel = 'OK',
                                                            content = c("BED files are 0-based and CNVxplorer works with 1-based data. Therefore, before the analysis, CNVxplorer adds 1 b.p to the start of the genomic interval.")
                 ),

                            width = 12,
                            collapsible = FALSE,
                            closable = FALSE,
                            fileInput("upload_bed_file", label = h5(strong("CNVs regions (.tsv file)")), 
                                      accept = c(".bed", '.tsv'))),
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
        tablerCard(title = 'Top 10 most frequent HP terms (gene annotation)',
                   plotOutput('overview_hp_terms'),
                   width = 6),
        tablerCard(title = 'Top 10 most frequent anatomical entities (gene annotation)',
                   plotOutput('overview_hp_anatomy'),
                   width = 6)),

        fluidRow(
          tablerCard(title = 'Indicate here phenotype terms of the patient',
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
                 uiOutput('n_hp_chosen'),
                 
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

                          tablerCard(title = 'Gene-disease associations' %>% helper(type = "inline",
                                                                                    style = "text-indent: 0.5em;",
                                                                                    title = "Similarity score",
                                                                                    size = "m",
                                                                                    buttonLabel = 'OK',
                                                                                    content = c(
                                                                                      "<b>Meaning of the '-' symbol</b>",
                                                                                      "The '-' symbol corresponds to p-values considered as non-significant or that could not be calculated due to the absence of phenotypic information. If the user has not entered any HP terms, this simbol represents a NA value",
                                                                                      "<b>Similarity score column is empty</b>",
                                                                                      "If you can not see any scores, you should first enter the patient's clinical symptoms in the panel above called 'Indicate here phenotype terms of the patient'.",
                                                                                      "<b>Two phenotypic similarity scores</b>",
                                                                                      "For each row, CNVxplorer displays the phenotypic score (p-value) corresponding to gene and disease. The reason for this, it is because you can have multiple diseases for the same gene and vice versa. Therefore, a gene and a disease that are associated may have a different phenotypic annotation.")),
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

                          tablerCard(title = 'HP terms associated with genes',
                                     DTOutput('hpo_assoc_genes'),
                                     width = 6),
                          tablerCard(title = 'HP terms associated with diseases',
                                     DTOutput('hpo_assoc_diseases'),
                                     width = 6)
                 )),
          
          tablerCard(title = 'Anatomical entities associated with HP terms',
                     highchartOutput('plot_anatomy'),
                     width = 12),
          tablerCard(title = 'DECIPHER CNVs - Phenotypic similarity',
                     DTOutput('decipher_similarity'),
                     overflow = TRUE,
                     width = 12)
          # tablerCard(title = 'Phenotypic similarity score',
          #            
          #            plotOutput('plot_similarity_genes'),
          #            width = 12,
          #            options = tagList( 
          #              
          #              prettyRadioButtons(
          #                inputId = "select_sim_gene_disease",
          #                label = '', 
          #                choices = list('Genes' = 'genes', 
          #                               'OMIM diseases' = 'diseases', 
          #                               'DECIPHER CNVs' = 'decipher',
          #                               'DECIPHER CNVs / Overlap' = 'decipher_overlap'),
          #                inline = TRUE, 
          #                status = "primary",
          #                fill = TRUE
          #              ))
          # )

        )
      ),
      tablerTabItem(
        tabName = "model",
        tablerCard(title = 'Mouse phenotypes associated with genes',
                   plotOutput('agg_model'),
                   width = 12),
        

          tablerCard(title = 'Orthologous genes associated with phenotypes in mouse (MGI)',
                     DTOutput('model_genes'),
                     width = 12)
        
      ),
      tablerTabItem(
        tabName = "docu",
      
          tablerCard(
            collapsible = FALSE,
            closable = FALSE,
            title = 'Documentation',
            fluidRow(
            column(12, align="center",
                   
            prettyRadioButtons(
              inputId = "select_doc_element",
              label = '', 
              choices = list('Overview' = 'overview',
                             'Tutorials' = 'tutorials',
                             'FAQs' = 'faqs',
                             'Input data' = 'versions',
                             'Installation' = 'installation',
                             'Browser compatibility' = 'browser', 
                             'About' = 'contact',
                             'Terms of use' = 'terms_of_use'),
              inline = TRUE, 
              status = "primary",
              fill = TRUE
            ))),
            uiOutput('doc_chosen'),
                     width = 12)
        
      ),
      tablerTabItem(
        tabName = "pubmed",
        fluidRow(
          tablerCard(
            title = "Frequency of genetic and phenotypic entities in titles and abstracts",
            
            collapsible = FALSE,
            closable = FALSE,
            plotOutput('pubtator_plot_disease') %>% withSpinner(type = 5),
            overflow = TRUE,
            width = 12,
            options = tagList(
              prettyRadioButtons(
                inputId = "select_del_dup_frequency",
                label = '',   
                choices =  split(c('deletions', 'duplications'), c('Deletions', 'Duplications')),
                inline = TRUE, 
                status = "primary",
                fill = TRUE
              )
            )
          ),
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
                        max = 50, value = 2),
            prettyRadioButtons(
              inputId = "entity_title_abstract",
              label = tags$b("Choose:"),
              choices = c('Title + abstract', "Title"),
              inline = FALSE, 
              status = "primary",
              fill = TRUE
            )
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
        
      )

    )
    
  )
)