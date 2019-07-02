


# Load libraries
library(shiny)
library(tidyverse)
library(tablerDash)
library(shinyEffects)
library(echarts4r)
library(shinyWidgets)
library(karyoploteR)
library(DT)
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
          tabName = "Home",
          icon = "home",
          "Home"
        ),
        tablerNavMenuItem(
          tabName = "Page1",
          icon = "box",
          "Overview"
        ),
        tablerNavMenuItem(
          tabName = "fa",
          icon = "box",
          "Functional analysis"
        ),
        tablerNavMenuItem(
          tabName = "model",
          icon = "box",
          "Model organism"
        ),
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
          tabName = "gene",
          icon = "plus",
          "Gene Panel"
        ),
        tablerNavMenuItem(
          tabName = "down_report",
          icon = "download",
          "Download report"
        ),
        tablerNavMenuItem(
          tabName = "Documentation",
          icon = icon("doc", lib = 'glyphicon'),
          "Documentation"
        )),
      id = "mymenu",
      src = "https://preview.tabler.io/demo/brand/tabler.svg",
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
    title = "CNVxplore",
    body = tablerDashBody(
      
      tablerTabItems(
        
        tablerTabItem(
          tabName = "Home",
          fluidRow(
          tablerCard(
            title = "Welcome to CNVxplore",
            
            'Description and more',
            width = 9,
            overflow = TRUE
          ),
          tablerCard(width = 3,
          tablerTimeline(
            tablerTimelineItem(
              title = "Overview",
              status = "green",
              date = ""
            ),
            tablerTimelineItem(
              title = "Functional analysis",
              status = NULL,
              date = ""
            ),
            tablerTimelineItem(
              title = "Model organism",
              status = NULL,
              date = ""
            ),
            tablerTimelineItem(
              title = "Regulatory regions",
              status = NULL,
              date = ""
            ),
            tablerTimelineItem(
              title = "Tissue-specificity",
              status = NULL,
              date = ""
            )
          ))),
          tablerCard(
            title = "Welcome to CNVxplore",
            
           DTOutput('score_references'),
            width = 12,
            overflow = FALSE
          )
        ),

        tablerTabItem(
          tabName = "Page1",
      
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
          # uiOutput('n_genes_pli'),
          # uiOutput('n_genes_rvis'),
          
          uiOutput('n_genes'),
          uiOutput('n_enhancer'),
          
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
          title = "Gene dataset",
 
          DTOutput("dgenes"),
          width = 12,
          overflow = TRUE
        ),
       
          column(width = 3,

          uiOutput('n_clinvar'),
          uiOutput('n_gwas'),
          uiOutput('n_lncrna'),
          uiOutput('n_hi')),
          column(width = 3,
                 
          uiOutput('n_dev'),
          uiOutput('n_tads'),
          uiOutput('n_omim'),
          uiOutput('n_pubmed'),
          uiOutput('n_pli')),
        column(width = 6,
        tablerCard(
          title = "Genome-wide percentile (pLI - RVIS scores)",
          
          plotlyOutput("dot_comparison"),
          width = 12,
          overflow = TRUE
        ))
      
    ),
    fluidRow(
      tablerCard(
        title = "Heatmap - Genome-wide percentile",
        plotlyOutput('heatmap'),
        width = 12,
        overflow = TRUE
        # options = tagList(
        #   pickerInput(
        #           inputId = "Id083",
        #           label = "Choose variables:",
        #           choices = list('pLI' = 'pLI'),
        #           multiple = TRUE
        #         )
        #   
        #   
        # )
      )
      
      
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
        title = "Gene ontology",
        width = 12,
        zoomable = FALSE,
        closable = FALSE,
        uiOutput('plot_df') ,
        options = tagList(
          switchInput(
            inputId = "enable_group_go",
            label = "Run?",
            value = FALSE,
            onStatus = "success",
            offStatus = "danger"
            
          ),
          switchInput(
            inputId = "user_df_plot",
            label = "Table?",
            value = FALSE,
            onStatus = "success",
            offStatus = "danger"
            
          ),
          pickerInput(
            inputId = "user_level",
            label = tags$b("Level:"), 
            choices = 1:52,
            width = 130,
            options = list(
              size = 10,
              `live-search` = TRUE)),
          prettyRadioButtons(
            inputId = "choose_group_go",
            label = tags$b("Select one option:"), 
            choices = c("biological process", "molecular function", "cellular component"),
            inline = TRUE, 
            status = "primary",
            fill = TRUE
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
        zoomable = FALSE,
        closable = FALSE,
        plotOutput('func_analysis') %>% withSpinner(type = 5),
        options = tagList(
          switchInput(
            inputId = "enable_func_analysis",
            label = "Run?",
            value = FALSE,
            onStatus = "success",
            offStatus = "danger"
            
          ),
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
          )
        )),
      tablerCard(
        title = "Pathway analysis",
        width = 12,
        zoomable = FALSE,
        closable = FALSE,
        plotOutput('func_pathways') %>% withSpinner(type = 5)
        ),
      tablerCard(
        title = "Gene-disease association",
        plotOutput('func_analysis_diseases') %>% withSpinner(type = 5),
        width = 12,
        zoomable = FALSE,
        closable = FALSE
      )
      
    ),
    tablerTabItem(
      tabName = "reg_region",
      fluidRow(
        tablerCard(title = 'Select a region:',
                   uiOutput('gen2e_2tissue'),
                   width = 3),
        uiOutput('n_enhancer_total'),
        uiOutput('n_enhancer_inside'),
        uiOutput('redund_n_enhancer'),
        tablerCard(title = 'List enhancers',
                   DTOutput('df_enhancer'),
                   width = 9)),
      tablerCard(title = 'RNA Expression (GTEx)',
                 plotlyOutput('tissue42_gtex'),
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
        tablerCard(title = 'Select a gene:',
                   uiOutput('n_pub2med'),
                   width = 3),
        tablerCard(title = 'Protein Expression (Human Protein Atlas)',
                   DTOutput('disease_pubmed'),
                   width = 9)),
      tablerCard(title = 'RNA Expression (GTEx)',
                 plotlyOutput('tissue2_555gtex'),
                 width = 12)
      
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
  tabName = "down_report",
  fluidRow(
    tablerCard(
      title = "Funnel overview",
      
      echarts4rOutput("funnel_genes"),
      width = 5,
      overflow = TRUE
    ),
    tablerCard(title = '',
               DTOutput('mode2l_genes'),
               width = 9))
  
)


      
    
    
    
    
    
    )
   
    )
  ),
  server = function(input, output) {
    
    gene_selected <- reactive({

      test4 <<- input$dgenes_rows_selected
      print(test4)
    })
    
    
    
    
    output$gene_tissue <- renderUI({
      
      input_data <- data_selected() %>% select(gene) %>% pull()
      
      pickerInput(
        inputId = "input_gene_tissue",
        # label = "Select gene:", 
        choices = input_data,
        options = list(
          size = 10,
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
    
    
    query_pubmed <- reactive({
      
      query_pubmed <- entrez_search(db="pubmed", term="22q11.2", retmax = 200 )
   

      
    })
    
    output$n_pubmed <- renderUI({
      
      tablerStatCard(
        value =   length(query_pubmed()[['ids']]),
        title = "Number of articles found in Pubmed",
        width = 12
      )
      
    })
    
    
    output$disease_pubmed <- renderDT({
      
      test21 <<- query_pubmed()
      query_tmp <- entrez_summary(db="pubmed", id= query_pubmed()[['ids']])
      
      title <- unname(map_chr(query_tmp, function(x) x[["title"]]))
      n_cites <- unname(map_chr(query_tmp, function(x) x[["pmcrefcount"]]))
      
      df_output <- tibble(title = title, n_cites = n_cites)

  
      datatable(df_output, rownames = FALSE, filter = 'top', 
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                  selection = 'single'
                  # columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
      
    })
    
    
    
    data_selected <- reactive({
      
      req(input$start_analysis > 0)
      
      data_raw <- hgcn_genes %>% filter(chrom == input$input_chrom)
      
      if(input$input_geno_karyo == 'Genomic coordinates') {
        
        data_raw <- data_raw  %>% mutate(keep = NA)
        for (i in 1:nrow(data_raw)) {
          data_raw$keep[i] <- c(data_raw$start_position[i], data_raw$end_position[i]) %overlaps% c(input$int_start, input$int_end)

        }

        data_raw <- data_raw %>% filter(keep == TRUE) %>% select(-keep)

      } else {
        
        data_raw  <- data_raw %>% filter(location == paste0(input$input_chrom, input$input_karyotype))

      }
      
      data_raw <- data_raw %>% 
        rename(band = location) %>% 
        mutate(coordinates = paste0(chrom,':', start_position,'-', end_position)) %>% 
        mutate(oe = paste0(oe_lof, ' (', oe_lof_lower, '-',oe_lof_upper, ')')) %>%
        select(-ensembl_gene_id, -transcript, -oe_lof, -oe_lof_lower, -oe_lof_upper, -vg)
      
      test1 <<- data_raw
    })
    
    output$dgenes <- renderDataTable({
      
      server <- TRUE
      data_input <- data_selected() %>% select(-start_position, -end_position) %>% select(-chrom)

      datatable(data_input, rownames = FALSE, filter = 'top', 
                # extensions = 'Responsive',
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                  selection = 'single',
                  columnDefs = list(list(className = 'dt-center', targets = '_all'))
                )) %>%
        formatStyle(
          'pLI',
          background = styleColorBar(c(0,1), '#ca7171'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })
    
    output$score_references <- renderDataTable({
      

      datatable(ref_scores, rownames = FALSE,
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap',
                  selection = 'single',
                  columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
    })
    
    output$distPlot <- renderPlot({
      if (input$enable_distPlot) hist(rnorm(input$obs))
    })
    
    output$tissue_gtex <- renderPlotly({
 
      
      filtered_gene <- input$input_gene_tissue
      p <- gtex %>%
        filter(gene == !!filtered_gene) %>%
        ggplot(aes(reorder(tissue, -value), value)) +
        geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
        theme_minimal() +
        xlab('Tissue') +
        ylab(paste('Median TPM')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position='none') +
        ggtitle(paste('Gene expression:', filtered_gene))
      
      ggplotly(p)
      
    })
    
    output$tissue_hpa <- renderDT({
      
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
        theme_minimal() +
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
        theme_minimal() +
        gghighlight(gene == name_gene_filtered)
      
      ggplotly(p)
    })
    
    
    
    
    
    
    output$data <- renderTable({
      mtcars[, c("mpg", input$variable), drop = FALSE]
    }, rownames = TRUE)
    
    output$plot_chrom <- renderPlot({
      
      input_chr <- paste0('chr', input$input_chrom)
      
      # data_input <- data_selected()
      # plotKaryotype()
      
      # ideoTrack <- IdeogramTrack(genome="hg19", chromosome= input$input_chr)
      # plotTracks(ideoTrack, from= input$int_start , to= input$int_end, showBandId=TRUE,
      #            cex.bands=0.5)
      # 
      ideoTrack <- IdeogramTrack(genome="hg19", chromosome= 'chr1')
      plotTracks(ideoTrack, from= 1000 , to= 10000000, showBandId=TRUE,
                 cex.bands=0.5)
        
      # plotKaryotype(chromosomes = input_chr, plot.type = 2) %>%
      # kpDataBackground(data.panel = 1)
      # kpAddBaseNumbers() %>%
      # kpRect(chr= input_chr, x0= data_input$start_position, x1= data_input$end_position, y0=0.2, y1=0.4)
      
      
      
    })
    
    output$choose_geno_karyo1 <- renderUI({
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        numericInput(
          inputId = "int_start",
          label = "Genomic interval - Start",
          value = 1000)

      } else {
        
        karyotype_filtered <- as.list(chromPlot::hg_cytoBandIdeo %>% filter(Chrom == input$input_chrom) %>% select(Name))

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
        
        name_region <- paste0('chr',input$input_chrom, ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(input$input_chrom, input$input_karyotype)
        
      }
      
      test12 <<- name_region
      
      tablerStatCard(
        value = nrow(data_selected()),
        title = HTML(paste0("Number of genes<br/>", name_region)),
        # trend = 19192,
        width = 12
      )
    })
    
    
    output$ref_user_genes <- renderUI({
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        name_region <- paste0('chr',input$input_chrom, ':', input$int_start, '-', input$int_end)
        
      } else {
        name_region <- paste0(input$input_chrom, input$input_karyotype)
        
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
      
      req(input$start_analysis > 0)
      
      
      if (input$input_geno_karyo == 'Genomic coordinates') {
        
        length_region <-  round((input$int_end - input$int_start + 1) / 1e6,2)
        
      } else {
        
        # name_region <- paste0(input$input_chrom, input$input_karyotype)
        length_region <-  round((input$int_end - input$int_start + 1) / 1e6,2)
      }
      
      tablerInfoCard(
        width = 12,
        value = paste(length_region, " Mb"),
        status = "primary",
        icon = "database",
        description =  'Length of the genomic region'
      )
      
      
    })
    
    
    output$n_lncrna <- renderUI({
      
      req(input$start_analysis > 0)
      
      data_tmp <- lncrna_coord %>% filter(chrom == input$input_chrom) %>%
        mutate(keep = 0)
      
      for (i in 1:nrow(data_tmp)) {
        data_tmp$keep[i] <- c(data_tmp$start[i], data_tmp$end[i]) %overlaps% c(input$int_start, input$int_end)
      }
      
      data_tmp <- data_tmp %>% filter(keep == 1) %>% select(id) %>% distinct() %>% pull(id)
  
        tablerStatCard(
        value =  length(data_tmp),
        title = "Number of lncRNA disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    prev_enhancer <- reactive({
      
      req(input$start_analysis > 0)
      
      data_tmp <- df_enhancers %>% filter(chrom == input$input_chrom) %>%
        mutate(keep = 0)
      
      for (i in 1:nrow(data_tmp)) {
        data_tmp$keep[i] <- c(data_tmp$start[i], data_tmp$end[i]) %overlaps% c(input$int_start, input$int_end)
      }
      data_tmp <- data_tmp %>% filter(keep == 1) %>% select(-keep)
    })
    
    output$n_enhancer_inside <- renderUI({
      
      list_genes_total <- prev_enhancer() %>% select(gene) %>% distinct() %>% pull(gene)
      data_tmp <- data_selected() %>% filter(gene %in% list_genes_total) %>% pull(gene)
      
      tablerStatCard(
        value =  paste(length(data_tmp), length(list_genes_total), sep = '/'),
        title = "Number of target genes inside CNV",
        width = 12
      )
    })
    
    
    
    output$n_enhancer <- renderUI({
      
      data_tmp <- prev_enhancer() %>% select(id) %>% distinct() %>% pull(id)
      

      tablerStatCard(
        value =  length(data_tmp),
        title = "Number of Enhancers",
        # trend = -10,
        width = 12
      )
    })
    
    
    redundancy_enhancers <- reactive({
      
      data_tmp <- prev_enhancer() %>% select(gene) %>% distinct() %>% pull(gene)
      df_output <- df_enhancers %>% filter(gene %in% data_tmp) %>% count(gene)
      test25 <<- df_output
    })
    
    output$redund_n_enhancer <- renderUI({
      
      test24 <<-  redundancy_enhancers() 

      data_tmp <- redundancy_enhancers() %>% filter(n == 1)


      tablerStatCard(
        value =  nrow(data_tmp),
        title = "Number of genes whose have one enhancer and it is disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    output$df_enhancer <- renderDT({
      
                datatable(prev_enhancer(), rownames = FALSE, filter = 'top', 
                          options = list(
                            pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                            selection = 'single'
                            # columnDefs = list(list(className = 'dt-center', targets = '_all'))
                          ))
      
      
    })
    
    output$n_tads <- renderUI({
      
      req(input$start_analysis > 0)
      
      
      n_tads <- check_tads(input$input_chrom, input$int_start, input$int_end )
      n_total <- nrow(data_selected())
      tablerStatCard(
        value =  length(n_tads),
        title = "Number of TADs disrupted",
        # trend = -10,
        width = 12
      )
    })
    
    output$n_dev <- renderUI({
      

      n_dev_yes <- data_selected() %>% filter(dev == 1) %>% nrow()
      n_total <- nrow(data_selected())
      tablerStatCard(
        value =  paste(n_dev_yes, n_total, sep = '/'),
        title = "Developmental disorder genes",
        # trend = -10,
        width = 12
      )
    })
    
    output$n_clinvar <- renderUI({
      
      
      n_clinvar_yes <- data_selected() %>% filter(clinvar == 1) %>% nrow()
      n_total <- nrow(data_selected())  
      
      tablerStatCard(
        value =  paste(n_clinvar_yes, n_total, sep = '/'),
        title =  'Genes located in ClinVar',
        # trend = -10,
        width = 12
      )
    })
    
    output$n_gwas <- renderUI({
      
      
      n_gwas_yes <- data_selected() %>% filter(gwas == 1) %>% nrow()
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
      
    gene_id <- data_selected() %>% select(gene) %>% pull()
    n_omim  <- morbidmap %>% filter(gene %in% gene_id) %>% nrow()
    n_total <- nrow(data_selected())  
    
     # morbidmap
      
      tablerStatCard(
        value =  paste(n_omim, n_total, sep = '/'),
        title = "OMIM diseases",
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
          value = 10000000)
      }
    })
    
    model_genes_phenotype <- reactive({
      
      
      go <- Ontology("mp")

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

      funnel <- data.frame(stage = c("19,192 genes", "182 genes", "60 genes"), value = c(1, 0.5, 0.25))
      
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
    
    
    output$heatmap <- renderPlotly({

      # a <- matrix(NA, nrow = 2, ncol = 100)
      # a[1,] <- hgcn_genes$pLI[1:100]
      # a[2,] <- hgcn_genes$rvis[1:100]
      # a[1,][as.numeric(which(is.na(a[1,])))] <- 0
      # a[2,][as.numeric(which(is.na(a[2,])))] <- 0
      # colnames(a) <- hgcn_genes$gene[1:100]
      # a <- pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE)
      # a
      n_genes <- nrow(data_selected())
      data_raw <- data_selected() %>% mutate(p_li = ntile(pLI, 100), p_rvis = ntile(rvis, 100),
                                        p_ncrvis = ntile(ncrvis, 100), p_ncgerp = ntile(ncgerp, 100))
      m <- matrix(NA, nrow = 4, ncol = n_genes)
      m[1,] <- data_raw$p_li[1:n_genes]
      m[2,] <- data_raw$p_rvis[1:n_genes]
      m[3,] <- data_raw$p_ncrvis[1:n_genes]
      m[4,] <- data_raw$p_ncgerp[1:n_genes]
      # m <- m[colSums(!is.na(m)) > 0]
      
      # heatmaply(as.matrix(m))
      
      plot_ly(
        x = data_selected() %>% select(gene) %>% pull(), y = c("RVIS", "pLI", 'ncRVIS', 'ncGERP'),
        z = m, 
        type = "heatmap"
        # width = 1200,
        # height = 500
      )
    })
    
    
    output$func_analysis <- renderPlot({
      
      if (input$enable_func_analysis == TRUE) {
      
        if (input$choose_go == 'biological process') {
          go_chosen <- 'BP'
        } else if (input$choose_go == 'molecular function') {
          go_chosen <- 'MF'
        } else {
          go_chosen <- 'CC'
        }
        
        filtered_genes <- data_selected() %>% select(entrez_id) %>% pull()  %>% as.character()
        univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character()

      go_analysis <- enrichGO(gene  = filtered_genes,
                      universe      = univ,
                      OrgDb         = org.Hs.eg.db,
                      ont           = go_chosen,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)
      
      validate(
        need(nrow(filtered_genes) != 0, "0 enriched terms found")
      )
      
      pathway_analysis <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
      
      a <-  go_analysis %>%
        as_tibble() %>%
        slice(1:10) %>%
        mutate(p.adjust = -log10(p.adjust)) %>%
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
        theme_minimal() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
      a

      
      }
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
      
      filtered_genes <- data_selected() %>% select(entrez_id) %>% pull()  %>% as.character()
      
      tryCatch(
        
        ggo <- groupGO(gene     = filtered_genes,
                       OrgDb    = org.Hs.eg.db,
                       ont      = go_chosen,
                       level    = n_level,
                       readable = TRUE),
        
        
        error= function(e) stop("Please, reduce the level assigned"))
      test20 <<- ggo
      
      
      validate(
        need(ggo %>% as_tibble() %>% select(Count) %>% sum() != 0, "0 terms found. Please reduce the level assigned")
      )
  
      ggo
      
    })
    
    
    output$plot_df <- renderUI({
      
      if (isTRUE(input$user_df_plot)) {
        
        DTOutput('df_go')
        
      } else {
        plotOutput('group_go')
        
      }

    })
    
    
    output$group_go  <- renderPlot({
      
      running_go() %>%
        as_tibble() %>%
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
        theme_minimal() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
      
      
    })
    
    
    output$df_go  <- renderDT({
      
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
    
    output$func_pathways  <- renderPlot({
      
      filtered_genes <- data_selected() %>% select(entrez_id) %>% pull()  %>% as.character()
      
      validate(
        need(nrow(filtered_genes) != 0, "0 enriched terms found")
      )

      pathway_analysis <- enrichPathway(gene= filtered_genes ,pvalueCutoff=0.05, readable=T)
      pathway_analysis %>%
        as_tibble() %>% 
        mutate(p.adjust = -log10(p.adjust)) %>%
        ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
        geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
        scale_fill_viridis_d() +
        coord_flip() +
        xlab('') +
        ylab('-log10(p-adjusted)') +
        ggtitle('Enriched pathways') +
        geom_hline(yintercept =-log10(as.numeric(input$sign_vline)), color = 'red', alpha = 0.6, linetype = 'dashed', size = 2) +
        geom_text(
          aes(label = GeneRatio, y = p.adjust + 0.05),
          position = position_stack(vjust = 0.5),
          vjust = 0
        ) +
        theme_minimal() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
    })
    
    output$func_analysis_diseases  <- renderPlot({
       
      # enrichNCG
      # enrichDGN
      # gseNCG
      # gseDGN
      
      test <<- data_selected()
      list_genes <- test  %>% slice(1:100) %>% select(entrez_id) %>% pull()  %>% as.character()
      
      # analysis_input <- enrichDGN(list_genes, minGSSize = 5)
      # barplot(analysis_input, showCategory=20)
      
      # enrich_dgn <- enrichDGN(list_genes)
      enrich_dgn <- enrichDO(list_genes)
      
      
      validate(
        need(nrow(enrich_dgn) != 0, "0 enriched terms found")
      )
      cnetplot(enrich_dgn, foldChange= hgcn_genes$entrez_id, readable = TRUE)
      
    })
    
    # output$func_analysis_mesh  <- renderPlot({
    #   
    #   # x <- enrichMeSH(de, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', category = 'C')
    #   dotplot(x)
    #   
    #   
    # })
    
  }
)
