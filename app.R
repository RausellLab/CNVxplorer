


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


human_chrom <- hg19 <- list('chr1' = 1, 'chr2' = 2,'chr3' = 3,'chr4' = 4,'chr5' = 5,'chr6' = 6,'chr7' = 7,'chr8' = 8,'chr9' = 9,'chr10' = 10,'chr11' = 11,'chr12' = 12,'chr13' = 13,
                            'chr14' = 14,'chr15' = 15,'chr16' = 16,'chr17' = 17,'chr18' = 18,'chr19' = 19, 'chr20' = 20, 'chr21' = 21,  'chr22' = 22,
                            'chrX' = 'X','chrY' = 'Y')

hg_cytoBandIdeo <- chromPlot::hg_cytoBandIdeo

# datas flowGl
vectors <- expand.grid(x = -3:3, y = -3:3)
mu <- 1
vectors$sx <- vectors$y
vectors$sy <- mu * (1 - vectors$x^2) * vectors$y - vectors$x
vectors$color <- log10(runif(nrow(vectors), 1, 10))

# calendar plot
dates <- seq.Date(as.Date("2018-01-01"), as.Date("2018-12-31"), by = "day")
values <- rnorm(length(dates), 20, 6)

year <- data.frame(date = dates, values = values)

# cards
flowCard <- tablerCard(
  title = "Karyotype",
  closable = FALSE,
  zoomable = FALSE,
  options = tagList(
    tablerAvatar(status = "lime", url = "https://john-coene.com/img/profile.png"),
    tablerTag(name = "build", addon = "passing", addonColor = "success")
  ),
  width = 12,
  plotOutput("flowGl")
)



profileCard <- tablerProfileCard(
  width = 12,
  title = "Peter Richards",
  subtitle = "Big belly rude boy, million
  dollar hustler. Unemployed.",
  background = "https://preview.tabler.io/demo/photos/ilnur-kalimullin-218996-500.jpg",
  src = "https://preview.tabler.io/demo/faces/male/16.jpg",
  tablerSocialLinks(
    tablerSocialLink(
      name = "facebook",
      href = "https://www.facebook.com",
      icon = "facebook"
    ),
    tablerSocialLink(
      name = "twitter",
      href = "https://www.twitter.com",
      icon = "twitter"
    )
  )
)



calendarCard <- tablerBlogCard(
  horizontal = TRUE,
  width = 12,
  echarts4rOutput("calendar")
)


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
          "Page1"
        ),
        tablerNavMenuItem(
          tabName = "Page2",
          icon = "box",
          "Page2"
        ),
        tablerNavMenuItem(
          tabName = "Documentation",
          icon = "box",
          "Documentation"
        )),
      id = "mymenu",
      src = "https://preview.tabler.io/demo/brand/tabler.svg",
      tablerDropdown(
        tablerDropdownItem(
          title = "Item 1 title",
          href = "http://google.com",
          url = "https://image.flaticon.com/icons/svg/1301/1301804.svg",
          status = "danger",
          date = "now",
          "This is the first dropdown item"
        ),
        tablerDropdownItem(
          url = "https://image.flaticon.com/icons/svg/1301/1301809.svg",
          status = "warning",
          "This is the second dropdown item",
          date = "yesterday"
        ),
        tablerDropdownItem(
          title = "Item 3 title",
          "This is the third dropdown item"
        )
      )
    ),
    footer = tablerDashFooter(
      tablerIcon(name = "maestro", lib = "payment"),
      tablerIcon(name = "mastercard", lib = "payment"),
      copyrights = "@David Granjon, 2019"
    ),
    title = "tablerDash",
    body = tablerDashBody(
      
      tablerTabItems(
        
        tablerTabItem(
          tabName = "Home",
          
          tablerCard(
            title = "Welcome to CNVxplore",
            
            'Description and more',
            width = 12,
            overflow = TRUE
          ),
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
         uiOutput('n_genes')
        ),
        column(6,
               plotOutput('flowGl')),
        column(
          width = 3,
          tablerCard(
            width = 12,
            tablerTimeline(
              tablerTimelineItem(
                title = "Gene dataset",
                status = "green",
                date = "now"
              ),
              tablerTimelineItem(
                title = "Genome-wide percentiles",
                status = NULL,
                date = "yesterday"
              ),
              tablerTimelineItem(
                title = "Functional analysis",
                status = NULL,
                date = "now"
              ),
              tablerTimelineItem(
                title = "Regulatory regions",
                status = NULL,
                date = "now"
              )
              )
            ),
          
          # uiOutput('n_genes_pli'),
          # uiOutput('n_genes_rvis'),
          

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
        fluidRow(
          uiOutput('p_pli'),
          uiOutput('p_rvis'),
          tablerStatCard(
            value =  '-',
            title = "Percentile pLI score",
            width = 2
          ),
          tablerStatCard(
            value =  '-',
            title = "Percentile pLI score",
            width = 2
          ),

          tablerStatCard(
            value =  '-',
            title = "Percentile pLI score",
            width = 2
          )),
        tablerCard(
          title = "Genome-wide percentile",
          
          plotOutput("perc_gene"),
          width = 6,
          overflow = TRUE
        ),
        tablerCard(
          title = "Genome-wide percentile (pLI - RVIS scores)",
          
          plotlyOutput("dot_comparison"),
          width = 6,
          overflow = TRUE
        )
      
    ),
    fluidRow(
      tablerCard(
        title = "Heatmap - Genome-wide percentile (pLI - RVIS scores)",
        plotlyOutput('heatmap'),
        width = 12,
        overflow = TRUE
      )
      
      
    ),
    tablerCard(
      title = "Functional analysis",
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
        prettyRadioButtons(
          inputId = "choose_go",
          label = "Select one option:", 
          choices = c("biological process", "molecular function", "cellular component"),
          inline = TRUE, 
          status = "primary",
          fill = TRUE
        )
      ))
        ))
   
    )
  ),
  server = function(input, output) {
    
    gene_selected <- reactive({

      test4 <<- input$dgenes_rows_selected
      print(test4)
    })
    
    
    
    data_selected <- reactive({
      
      data_raw <- hgcn_genes %>% filter(chrom == 1)
      
      if(input$input_geno_karyo == 'Genomic coordinates') {
        
        data_raw <- data_raw  %>% mutate(keep = NA)
        for (i in 1:nrow(data_raw)) {
          data_raw$keep[i] <- c(data_raw$start_position[i], data_raw$end_position[i]) %overlaps% c(input$int_start, input$int_end)

        }
        test1 <<- data_raw
        
        data_raw <- data_raw %>% filter(keep == TRUE) %>% select(-keep)

      } else {
        
        data_raw  <- data_raw %>% filter(location == paste0(input$input_chrom, input$input_karyotype))

      }
      
      data_raw <- data_raw %>% 
        rename(band = location) %>% 
        mutate(coordinates = paste0(chrom,':', start_position,'-', end_position)) %>% 
        mutate(oe = paste0(oe_lof, ' (', oe_lof_lower, '-',oe_lof_upper, ')')) %>%
        select(-ensembl_gene_id, -chrom, -entrez_id, -transcript, -oe_lof, -oe_lof_lower, -oe_lof_upper, -vg)

    })
    
    output$dgenes <- renderDataTable({
      
      server <- TRUE
      data_input <- data_selected() %>% select(-start_position, -end_position)
      datatable(data_input, rownames = FALSE, filter = 'top',
                options = list(
                  pageLength = 5, autoWidth = TRUE, style = 'bootstrap',
                  selection = 'single',
                  columnDefs = list(list(className = 'dt-center', targets = '_all'))
                ))
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
    
    
    output$perc_gene <- renderPlot({
     
      red_line <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(pLI) %>% pull()
      symbol_chosen <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(gene) %>% pull()
      
      hgcn_genes %>% 
        ggplot(aes(pLI)) + 
        geom_histogram() +
        theme_minimal() +
        scale_fill_viridis_d() +
        geom_vline(xintercept = red_line, color = 'red', alpha = 0.6, type = 'dashed') +
        ggtitle(paste('Histogram pLI scores - Gene:', symbol_chosen))
      
    })
    
    output$dot_comparison <- renderPlotly({
      
      name_gene_filtered <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(gene) %>% pull()
      test6 <<- name_gene_filtered
      test7 <<- data_selected()
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
    
    output$flowGl <- renderPlot({
      input_chr <- paste0('chr', input$input_chrom)
      data_input <- data_selected()
      # plotKaryotype(chromosomes = input_chr, plot.type = 2) %>%
      # kpDataBackground(data.panel = 1) %>%
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
        
        karyotype_filtered <- as.list(chromPlot::hg_cytoBandIdeo %>% filter(Chrom == input$chrom) %>% select(Name))
        
        
        selectizeInput(inputId = 'input_karyotype', label = 'Karyotype', choices = karyotype_filtered,
                       selected = NULL, multiple = FALSE,
                       options = NULL)

      }

    })
    
    output$n_genes <- renderUI({
      
      tablerStatCard(
        value = nrow(data_selected()),
        title = "Number of genes",
        trend = 19192,
        width = 12
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
    
    output$n_genes <- renderUI({
      
      tablerStatCard(
        value = nrow(data_selected()),
        title = "Number of genes",
        trend = -10,
        width = 12
      )
    })
    
    
    output$choose_geno_karyo1 <- renderUI({

        numericInput(
          inputId = "int_start",
          label = "Genomic interval - Start",
          value = 1000)

      
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
          value = 60000000)
      }
    })
    
    
    output$info <- renderUI({
      tablerInfoCard(
        width = 12,
        value = paste0(input$totalStorage, "GB"),
        status = "success",
        icon = "database",
        description = "Total Storage Capacity"
      )
    })
    
    
    output$heatmap <- renderPlotly({

      a <- matrix(NA, nrow = 2, ncol = 100)
      a[1,] <- hgcn_genes$pLI[1:100]
      a[2,] <- hgcn_genes$rvis[1:100]
      # a[1,][as.numeric(which(is.na(a[1,])))] <- 0
      # a[2,][as.numeric(which(is.na(a[2,])))] <- 0
      # colnames(a) <- hgcn_genes$gene[1:100]
      # a <- pheatmap(a, cluster_rows = FALSE, cluster_cols = FALSE)
      # a
      
      data_raw <- hgcn_genes %>% mutate(p_li = ntile(pLI, 100), p_rvis = ntile(rvis, 100),
                                        p_ncrvis = ntile(ncrvis, 100), p_ncgerp = ntile(ncgerp, 100))
      m <- matrix(NA, nrow = 4, ncol = 500)
      m[1,] <- data_raw$p_li[1:500]
      m[2,] <- data_raw$p_rvis[1:500]
      m[3,] <- data_raw$p_ncrvis[1:500]
      m[4,] <- data_raw$p_ncgerp[1:500]
      
      p <- plot_ly(
        x = hgcn_genes$gene[1:500], y = c("RVIS", "pLI", 'ncRVIS', 'ncGERP'),
        z = m, type = "heatmap"
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
      test <- hgcn_genes %>% slice(1:100) %>% select(entrez_id) %>% pull()  %>% as.character()
      univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character

      go_analysis <- enrichGO(gene  = test,
                      universe      = univ,
                      OrgDb         = org.Hs.eg.db,
                      ont           = go_chosen,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
      
      pathway_analysis <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
      
      
      a <-  go_analysis %>% ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
        geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
        scale_fill_viridis_d() +
        # scale_y_log10() +
        coord_flip() +
        xlab('') +
        ylab('p-adjusted') +
        ggtitle('GO analysis') +
        geom_vline(xintercept = 0.05) +
        theme_minimal() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
      b <-  pathway_analysis %>% ggplot(aes(reorder(Description, p.adjust), p.adjust)) +
        geom_col(aes(fill = Description), color = 'black', show.legend = FALSE) +
        scale_fill_viridis_d() +
        # scale_y_log10() +
        coord_flip() +
        xlab('') +
        ylab('p-adjusted') +
        ggtitle('Pathway analysis') +
        geom_vline(xintercept = 0.05) +
        theme_minimal() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      
      a + b
      }
    })
    
    output$path_analysis <- renderPlot({
      
      test <- hgcn_genes %>% slice(1:100) %>% select(entrez_id) %>% pull()  %>% as.character()
      univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character
      
      ego <- enrichGO(gene          = test,
                      universe      = univ,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
      
      barplot(ego, showCategory=30)
      
      
      
    })
    
  }
)
