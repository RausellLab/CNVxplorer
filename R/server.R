function(input, output, session) {
  
  
  observe_helpers() # do not remove -> helper dialogs 
  
  observeEvent(input$start_analysis, {
    
    shinyjs::reset("dgenes_rows_selected")
    shinyjs::reset("dgenes_rows_all")
    shinyjs::reset("chosen_hp")
    shinyjs::reset("data_selected_prev")
    shinyjs::reset("start_analysis")
    #
    shinyjs::reset('enhancers_on_off')
    shinyjs::reset('data_selected_enhancers')
    
    shinyjs::reset('tads_on_off')
    shinyjs::reset('data_selected_tads')
    
    shinyjs::reset('tfs_on_off')
    shinyjs::reset('data_selected_tfs')
    
    shinyjs::reset('mirnas_on_off')
    shinyjs::reset('data_selected_mirnas')
    
    shinyjs::reset('lncrnas_on_off')
    shinyjs::reset('data_selected_lncrnas')
    #
    shinyjs::reset('filter_by_gene_ppi')
    
    shinyjs::reset('input_gene_tissue')
    shinyjs::reset('input_tissue')
    #
    shinyjs::reset('counter_header')
    shinyjs::reset('select_reg_region')
    #
    shinyjs::reset('enable_net')
    shinyjs::reset('choose_net_source')
    #
    shinyjs::reset("input_inheritance")
    shinyjs::reset("hpo_filter_genes_rows_selected")
    shinyjs::reset("hpo_filter_diseases_rows_selected")
    #
    shinyjs::reset("enable_func_analysis")
    shinyjs::reset("enable_path_analysis")
    shinyjs::reset("enable_do_analysis")
    shinyjs::reset("enable_group_go") 
    shinyjs::reset("enable_tsea")
    #
    shinyjs::reset('only_omim')
    shinyjs::reset('select_del_dup')
    #
    shinyjs::reset('chosen_hp')
    
  })
  
  
  observeEvent(input$enhancers_on_off, {
    
    shinyjs::reset("enable_func_analysis")
    shinyjs::reset("enable_path_analysis")
    shinyjs::reset("enable_do_analysis")
    shinyjs::reset("enable_group_go")
    shinyjs::reset("enable_tsea")
    shinyjs::reset('chosen_hp')
    
  })
  
  observeEvent(input$tads_on_off, {
    
    shinyjs::reset("enable_func_analysis")
    shinyjs::reset("enable_path_analysis")
    shinyjs::reset("enable_do_analysis")
    shinyjs::reset("enable_group_go")
    shinyjs::reset("enable_tsea")
    shinyjs::reset('chosen_hp')
    
  })
  
  observeEvent(input$mirnas_on_off, {
    
    shinyjs::reset("enable_func_analysis")
    shinyjs::reset("enable_path_analysis")
    shinyjs::reset("enable_do_analysis")
    shinyjs::reset("enable_group_go")
    shinyjs::reset("enable_tsea")
    shinyjs::reset('chosen_hp')
    
  })
  
  observeEvent(input$tfs_on_off, {
    
    shinyjs::reset("enable_func_analysis")
    shinyjs::reset("enable_path_analysis")
    shinyjs::reset("enable_do_analysis")
    shinyjs::reset("enable_group_go")
    shinyjs::reset("enable_tsea")
    shinyjs::reset('chosen_hp')
    
    
  })
  
  observeEvent(input$input_geno_karyo, {
    
    shinyjs::reset('int_start')
    shinyjs::reset('int_end')
    shinyjs::reset('input_chrom')
    shinyjs::reset('file_cnv') 
    shinyjs::reset('cnv_file')
  })
  
  
  
  observeEvent(input$reset_pheno_analysis, {
    
    
    shinyjs::reset('chosen_hp')
  })
  
  observeEvent(input$input_geno_karyo, {
    
    
    shinyjs::reset('file_cnv')
  })
  
# Check input amenabar
  
  observe({

    test_negative <- nrow(coord_user() %>% filter(start < 0 | end < 0))
    # negative numbers
    if (test_negative > 0) {

      shinyalert("Error!", "There is a negative coordinate", type = "error")

      req(coord_user() %>% pull(start) > 0)


    # out of chromosome coordinate
    } else if (nrow(coord_user() %>%
               left_join(coord_chrom_hg19, by = 'chrom') %>%
               mutate(is_bigger = end - length) %>%
               filter(is_bigger > 0)
               ) > 0) {
      shinyalert("Error!", "You have selected a genomic coordinate out of the chromosome", type = "error")
    } else if (nrow(coord_user() %>% mutate(t_length = end - start + 1) %>%
                    filter(t_length > 1.5e7))) {

      shinyalert("Error!", "One of the genomic intervals exceed the maximum length (10 millions b.p)", type = "error")

    }
  })

  map_blacklist <- reactive({

    input_test <- coord_user() %>% mutate(start = start - 1) 
    
    tmp_df <- blacklist_encode %>% 
      mutate(length_cnv = end - start + 1) %>% 
      bed_coverage(input_test) %>%
      filter(.cov > 0) %>%
      mutate(.frac = paste0(round(.frac * 100, 2), '%')) %>%
      mutate(tag = paste0(class,' - ', .frac, ' (', chrom,':', start, '-', end, ')'))
    
    
    tmp_df

  })
  
  
  output$check_blacklist <- renderUI({
    
    tmp_df <- map_blacklist() 
    
    
    req(nrow(tmp_df) > 0)

    tablerInfoCard(
      width = 12,
      value =  paste(nrow(tmp_df), 'problematic regions'),
      status = "danger",
      icon = "database",
      description =  HTML(tmp_df %>% pull(tag) %>% paste(collapse = ' <br/> '))
    )
    
    
  })

  observeEvent(input$take_intersect, {
    
    
    df_tmp <- intersection_running() %>%
      select(id, start, end) %>%
      slice(input$df_intersection_rows_selected)
    
    new_start <- df_tmp %>% pull(start)
    new_end <- df_tmp %>% pull(end)
    
    
    updateNumericInput(session, "int_start", value = new_start)
    updateNumericInput(session, "int_end", value = new_end)
    
  })
  
  
  coord_user <- eventReactive(input$start_analysis, {

    coord_start <- as.numeric(str_remove_all(input$int_start, ','))
    coord_end <-  as.numeric(str_remove_all(input$int_end, ','))
    coord_chrom <- input$input_chrom


    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      
      tbl_output <- tibble('chrom' = coord_chrom, 'start' = coord_start,
                           'end' = coord_end)
      
      
    } else if (input$input_geno_karyo == 'G banding') {
      
      df_tmp <- coord_cytobands %>%
        filter(Chrom == coord_chrom) %>%
        filter(Name == input$input_karyotype)
      
      coord_start <- df_tmp %>% select(Start) %>% pull()
      coord_end <-  df_tmp %>% select(End) %>% pull()
      coord_chrom <- input$input_chrom
      
      coord_start <- as.numeric(coord_start)
      coord_end <- as.numeric(coord_end)
      
      tbl_output <- tibble('chrom' = coord_chrom, 'start' = coord_start,
                           'end' = coord_end)
      
    } else {

      tbl_output <- cnv_file_to_analyze()

    }
    

    tmp_check_1 <- tbl_output %>% 
      filter(start < 0 | end < 0)
    
    tmp_check_2 <- tbl_output %>% 
      mutate(t_length = end - start + 1) %>%
      filter(t_length >= 1.5e7)
    
    tmp_check_3 <- tbl_output %>%
      mutate(check_lower_end = end - start) %>%
      filter(check_lower_end < 0) 
    
    tmp_check_4 <- tbl_output %>%
      left_join(coord_chrom_hg19, by = 'chrom') %>%
      mutate(is_bigger = end - length) %>%
      filter(is_bigger > 0)
    
    
    if (nrow(tmp_check_1) > 0) {

      shinyalert("Error!", "You have selected a region with negative coordinates", type = "error")
      req(tmp_check_1 == 0)

    } else if (nrow(tmp_check_2) > 0) {
      
      shinyalert("Error!", "One of the genomic intervals exceed the maximum length (15 millions b.p)", type = "error")
      req(tmp_check_2 == 0)

    } else if (nrow(tmp_check_3) > 0 ) {
       
      
      shinyalert('Error!', "The end of the genomic interval is lower than the start", type = "error")
      req(tmp_check_3 == 0)
      
     } else if (nrow(tmp_check_4) > 0 ) {
 
      shinyalert("Error!", "You have selected a genomic coordinate out of the chromosome", type = "error")
      req(tmp_check_4 == 0)
       
      }
    
    tbl_output
  })
  
  
  
  
  output$gtex_gene <- renderUI({
    
    # if (is.null(input$dgenes_rows_all)) {
    #   df_genes <- data_selected()
    # } else {
    #   df_genes <- data_selected()[input$dgenes_rows_all,]
    # }
    
    input_data <- data_selected() %>% select(gene) %>% pull()
    
    pickerInput(
      inputId = "input_gtex_gene",
      # label = "Select gene:", 
      choices = input_data,
      options = list(
        size = 5,
        `live-search` = TRUE))
    
    
    
  })
  
  output$gtex_tissue <- renderUI({
    
    # if (is.null(input$dgenes_rows_all)) {
    #   df_genes <- data_selected()
    # } else {
    #   df_genes <- data_selected()[input$dgenes_rows_all,]
    # }
    
    input_data <- data_selected() %>% select(gene) %>% pull()
    
    input_tissue <- gtex %>% 
      filter(gene %in% input_data) %>% 
      select(tissue) %>% 
      distinct() %>% 
      pull(tissue)
    
    pickerInput(
      inputId = "input_gtex_tissue",
      # label = "Select gene:", 
      choices = input_tissue,
      options = list(
        size = 5,
        `live-search` = TRUE))

  })
  
  
  output$gene_tissue <- renderUI({
 
    input_data <- running_hpa() %>% select(gene) %>% distinct() %>% pull()
    
    pickerInput(
      inputId = "input_gene_tissue",
      choices = input_data,
      options = list(
        size = 5,
        `live-search` = TRUE))
    
    
    
  })
  
  output$select_tissue <- renderUI({
    
    # req(input$input_gene_tissue)
    
    choices_vector <- hpa %>% select(tissue) %>% distinct() %>% pull() %>% as.character()
    
    pickerInput(
      inputId = "input_tissue",
      # label = "Select gene:", 
      choices = choices_vector,
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
    
    
    df_output <-  cnv_df %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      group_by(chrom, start, end) %>%
      filter(.overlap == max(.overlap)) %>%
      ungroup() %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      distinct()
    

    df_output
    
  })

  check_hp_genes <- reactive({

    validate(
      need(!is.null(input$chosen_hp), "Please, select phenotype terms.")
    )
    
    
    hp_chosen <- input$chosen_hp
    genes_chosen <- df_genes %>% select(entrez_id) %>% pull()
    
    
    df_tmp <- hpo_genes %>%
      filter(hp %in% hp_chosen) %>%
      filter(entrez_id %in% genes_chosen)
    
    
    df_tmp2 <- df_tmp %>%
      select(hp) %>%
      distinct()

    
    df_tmp <- df_tmp %>% 
      left_join(df_tmp2, by = 'hp')
    
    df_tmp
    
    
  })
  
  
  suggest_hpo_terms <- reactive({
    
    validate(
      need(length(input$chosen_hp) > 0, "Please, select at least one HPO term.")
    )
    
    hp_chosen <- input$chosen_hp # input$chosen_hp
    hp_descendants <- get_descendants(hpo_dbs, hp_chosen)
    hp_descendants <- hp_descendants[!hp_descendants %in% hp_chosen] 
    
    
    if (length(hp_descendants) > 5) {
      
      hp_descendants <- ontologySimilarity::descendants_IC(hpo_dbs) %>% as_tibble(rownames = 'term') %>%
        filter(term %in% hp_descendants) %>% 
        arrange(desc(value)) %>% slice(1:5) %>% select(-value) %>% pull(term)
      
    }
    
    
    hp_descendants <- hp_descendants %>% paste(collapse = '|')
  })
  
  output$suggest_df <- renderDT({
    
    validate(
      need(nchar(suggest_hpo_terms()) > 0, "Select more HPO terms.")
    )
    
    
    tmp_df <- vector_total_terms %>% 
      as_tibble() %>% 
      filter(str_detect(term, suggest_hpo_terms() )) %>% 
      rename(Term = value) %>%
      select(-term_desc) %>%
      mutate(term = paste0("<a href='", paste0('https://hpo.jax.org/app/browse/term/', term),"' target='_blank'>", term,"</a>"))
    
    validate(
      need(nrow(tmp_df) > 0, "Select more HPO terms.")
    )
    
    datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('HPO term', 'Description'), options = list(dom = 't'), class = 'cell-border stripe')

  })
  
  
  get_entities <- reactive({
    
    req(input$text_recognition)

    b <- find_scigraph(text = input$text_recognition)
    
    b$content %>% 
      filter(category == 'phenotype', str_detect(identifier, 'HP:')) %>%
      select(-category)
  })
  
  output$entities_df <- renderDT({
    
    validate(
      need(nrow(get_entities()) > 0, 'No phenotypic terms found in the text.')
    )
    
    datatable(get_entities(), 
              rownames = FALSE, 
              options = list(dom = 't'), 
              class = 'cell-border stripe')
    
  })
  
  output$check_genes_hp <- renderUI({
    
    # req(input$start_analysis > 0)
    
    df_genes <- data_selected()
    
    
    
    tablerStatCard(
      value =  nrow(check_hp_genes()),
      title = "Number of genes associated with the phenotype(s)",
      width = 12
    )
    
  })


  
  output$df_check_hp_genes <- renderDT({
    
    # req(input$start_analysis > 0)
    

    validate(
      need(nrow(check_hp_genes()) != 0, "Not genes found.")
    )
    
    datatable(check_hp_genes(), options = list(
      pageLength = 5))
  })

  
  output$n_hp_chosen <- renderUI({
    
    # req(input$start_analysis > 0)
    
    tablerInfoCard(
      width = 12,
      value =  paste(length(input$chosen_hp), 'HPO terms'),
      status = "primary",
      icon = "database",
      description =  'We recommend more than three HP terms'
      
    )
    
    # tablerStatCard(
    #   value =  length(input$chosen_hp),
    #   title = "Clinical features selected",
    #   # trend = -10,
    #   width = 12
    # )
    
  })
  
  output$overview_hp_terms <- renderPlot({

    hpo_filter() %>% 
      count(hp) %>% 
      left_join(hpo_genes %>% select(hp, desc),by = 'hp') %>% 
      distinct() %>%
      mutate(desc = paste0(desc, ' (', hp, ')')) %>%
      arrange(desc(n)) %>%
      slice(1:10) %>%
      ggplot(aes(reorder(desc, n), n)) +
      geom_col(aes(fill = n), color = 'black', show.legend = FALSE) +
      scale_fill_viridis_c() +
      coord_flip() +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme_minimal() +
      labs(y = 'Frequency', x = NULL)

  })
  
  output$overview_hp_anatomy <- renderPlot({
    
    genes_selected <- hpo_filter() %>% select(gene) %>% distinct() %>% pull()
    hpo_from_gene <- hpo_genes %>% filter(gene %in% genes_selected)  %>% pull(hp)


   unlist(map(hpo_from_gene, function(x) get_ancestors(hpo_dbs, x))) %>% 
   enframe()  %>% 
   select(value) %>%
   left_join(anato_df, by = c('value' = 'name')) %>%
   na.omit() %>%
   rename(anatomy_entity = value.y) %>%
      count(anatomy_entity) %>%
      slice(1:10) %>%
      ggplot(aes(reorder(anatomy_entity, n), n)) +
      geom_col(aes(fill = n), color = 'black', show.legend = FALSE) +
      # scale_fill_continuous() +
      scale_fill_viridis_c() +
      coord_flip() +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme_minimal() +
      labs(y = 'Nº HP terms associated', x = NULL)

  })
  
  
  
  output$n_cnv_patho <- renderUI({
    
    # req(input$start_analysis > 0)
    
    data_tmp <- check_cnv_df() %>% 
      filter(source == 'decipher' & pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>% 
      nrow()
    
    data_tmp2 <- running_clinvar_yes_cnv() %>% nrow()
    
    tablerStatCard(
      value =  data_tmp + data_tmp2,
      title = "Pathogenic CNVs (DECIPHER & ClinVar)",
      width = 12
    )
  })
  
  output$ui_select_decipher_clinvar <- renderUI({
    
    req(running_clinvar_yes_cnv())
    req(check_cnv_df())
    
    data_tmp <- check_cnv_df() %>% 
      filter(source == 'decipher' & pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>% 
      nrow()
    
    data_tmp2 <- running_clinvar_yes_cnv() %>% nrow()
    
    
    vector_tmp <- c(paste('ClinVar',paste0('(', data_tmp2, ')')),
                    paste('DECIPHER', paste0('(', data_tmp, ')')))

    vector_n_dbs <-    rev(split(c('ClinVar', 'DECIPHER'), vector_tmp))
    
    prettyRadioButtons(
      inputId = "select_decipher_clinvar",
      label = '', 
      choices =  vector_n_dbs,
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
    
  })
  
  output$n_cnv_nopatho <- renderUI({
    
    # req(input$start_analysis > 0)
    
    data_tmp <- check_cnv_df() %>% filter(source != 'decipher') %>% nrow()
    
    tablerStatCard(
      value =  data_tmp,
      title = "Nonpathogenic CNVs (gnomAD, DGV, DECIPHER)",
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  data_tmp,
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of nonpathogenic CNVs",
    #   width = 12
    # )
    
    
  })
  
  
  query_pubmed_del <- reactive({
    
    
    if ((input$input_geno_karyo == 'Multiple coordinates (NGS)' | input$input_geno_karyo == 'Genomic coordinates')) {

      tmp_query <-  coord_cytobands %>%
        rename(chrom = Chrom, start = Start, end = End) %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        mutate(chrom_name = paste0(chrom, Name)) %>%
        select(chrom, Name, chrom_name) %>%
        mutate(result = paste('(','(', 'chromosome', chrom,'AND' , Name, ')' ,'OR', chrom_name, ')')) %>%
        pull(result) %>%
        paste(collapse = ' OR ') %>%
        paste('AND (deletion OR microdeletion) AND homo sapiens')
        
      
      query_region <- tmp_query
      

    } else {
      
      chrom_tmp <- paste('chromosome',  coord_user() %>%  pull(chrom))
      band_tmp <-  input$input_karyotype
      band2_tmp <- paste0( coord_user() %>%  pull(chrom), band_tmp)

      
      query_region <- paste('(','(', chrom_tmp,'AND', band_tmp,')','OR', band2_tmp,')', 'AND deletion AND homo sapiens')
      
    }
    
    query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 600 )

  })

  
  output$n_pubmed_del <- renderUI({

    tablerStatCard(
      value =   length(query_pubmed_del()[['ids']]),
      title = "Articles found in Pubmed associated with deletions",
      width = 12
    )
    
  })
  
  
  output$ui_only_omim <- renderUI({
    
    
    prettyRadioButtons(
      inputId = "only_omim",
      label = tags$b("Filter articles associated with OMIM entries:"), 
      choices = c("No", "Yes"),
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
    
    
  })

  omim_assoc <- reactive({
    
    
    if (input$select_del_dup == 'deletions') {
      
      ids_query <- query_pubmed_del()[['ids']]
    } else {
      ids_query <- query_pubmed_dup()[['ids']]
    }

    
    query_link <- entrez_link(db= 'omim', id= ids_query, dbfrom="pubmed", by_id = TRUE)

    tmp_df <- tibble(pubmed_id = ids_query, omim_assoc = NA)
    
    for (i in 1:length(query_link)) {
      
      tmp_value <-  query_link[[i]]$links$pubmed_omim_calculated
      tmp_value <- tmp_value %>% paste(collapse = ' ')
      
      if (is.null(tmp_value)) {
        next
      } else {
        
        tmp_df$omim_assoc[i] <- tmp_value
        
      }
      
    }
    
    tmp_df <- tmp_df %>% separate_rows(omim_assoc, sep = ' ')
    

    validate(
      need(nrow(tmp_df  %>% filter(omim_assoc != '')) != 0, 'No articles associated with OMIM entries.')
    )
    
    
    tmp_df
  })

  query_pubmed_dup <- reactive({

    if ((input$input_geno_karyo == 'Multiple coordinates (NGS)' | input$input_geno_karyo == 'Genomic coordinates')) {
      
      tmp_query <-  coord_cytobands %>%
        rename(chrom = Chrom, start = Start, end = End) %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        mutate(chrom_name = paste0(chrom, Name)) %>%
        select(chrom, Name, chrom_name) %>%
        mutate(result = paste('(','(', 'chromosome', chrom,'AND' , Name, ')' ,'OR', chrom_name, ')')) %>%
        pull(result) %>%
        paste(collapse = ' OR ') %>%
        paste('AND (duplication OR microduplication) AND homo sapiens')

      

      # if (length(tmp_query) > 1) tmp_query %>% coll 
      query_region <- tmp_query
      
      
    } else {
      
      chrom_tmp <- paste('chromosome',  coord_user() %>%  pull(chrom))
      band_tmp <-  input$input_karyotype
      band2_tmp <- paste0( coord_user() %>%  pull(chrom), band_tmp)
      
      query_region <- paste('(','(', chrom_tmp,'AND', band_tmp,')','OR', band2_tmp,')', 'AND duplication AND homo sapiens')
      
    }

    query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 600 )
    


  })
  
  output$n_pubmed_dup <- renderUI({
    
    
    tablerStatCard(
      value =   length(query_pubmed_dup()[['ids']]),
      title = "Articles found in Pubmed associated with duplications",
      width = 12
    )

  })
  
  
  output$n_mortality <- renderUI({
    
    if (nrow(model_genes_phenotype()) > 0) {
    tmp_n <- model_genes_phenotype() %>%  filter(description == 'mortality/aging') %>% nrow()
    } else {
    tmp_n <- 0
    }

    tablerStatCard(
      value =  tmp_n,
      title = "Genes associated with mortality/aging phenotype",
      width = 12
    )
    
    
  })
  
  output$n_embryo <- renderUI({
    
    if (nrow(model_genes_phenotype()) > 0) {
      tmp_n <- model_genes_phenotype() %>%
        filter(description == 'embryo phenotype') %>%
        nrow()
    } else {
      tmp_n <- 0
    }
    


    tablerStatCard(
      value =   tmp_n,
      title = "Genes associated with embryonic phenotype",
      width = 12
    )
    
  })

  
  
  running_pubmed_del <- reactive({
    
    req(coord_user())
    
    tryCatch(
        query_tmp <- entrez_summary(db="pubmed", id= query_pubmed_del()[['ids']]),
      error= function(e) stop("Server error, please try later.")
      )

    
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
    
    df_output <- df_output %>% 
      mutate(n_cites = as.integer(ifelse(n_cites == '', 0, n_cites)))
    
    df_output
    
    
  })
  
  running_pubmed_dup <- reactive({
    
    req(query_pubmed_dup())
    
    if (length(query_pubmed_dup()[['ids']]) == 0) {
      
      df_output <- tibble()
      
    } else {
    
    query_tmp <- entrez_summary(db="pubmed", id= query_pubmed_dup()[['ids']])
    
    if (length(query_tmp[[1]]) > 1) {
    
    title <- unname(map_chr(query_tmp, function(x) x[["title"]]))
    n_cites <- unname(map_chr(query_tmp, function(x) x[["pmcrefcount"]]))
    first_author <- unname(map_chr(query_tmp, function(x) x[["sortfirstauthor"]]))
    last_author <- unname(map_chr(query_tmp, function(x) x[["lastauthor"]]))
    date_release <- unname(map_chr(query_tmp, function(x) x[["pubdate"]]))
    pmid_id <- unname(map_chr(query_tmp, function(x) x[["uid"]]))
    journal_id <- unname(map_chr(query_tmp, function(x) x[["fulljournalname"]]))
    
    } else {
      
      title <- query_tmp[["title"]]
      n_cites <- query_tmp[["pmcrefcount"]]
      first_author <- query_tmp[["sortfirstauthor"]]
      last_author <- query_tmp[["lastauthor"]]
      date_release <- query_tmp[["pubdate"]]
      pmid_id <- query_tmp[["uid"]]
      journal_id <- query_tmp[["fulljournalname"]]
      
    }
    
    df_output <- tibble(title = title, first_author = first_author, last_author = last_author, 
                        n_cites = n_cites, journal = journal_id, date_release = date_release, 
                        pmid = pmid_id)
    
    df_output <- df_output %>% mutate(n_cites = as.integer(ifelse(n_cites == '', 0, n_cites)))
    
    }
    df_output
    
    
  })
  
  output$del_dup_pubmed <- renderDT({
    
    req(running_pubmed_del())
    req(running_pubmed_dup())
    req(input$select_del_dup)
    req(input$only_omim)
    
    
    if (input$select_del_dup == 'deletions') {
      
      tmp_df <- running_pubmed_del()
      
    } else {
      
      validate(
        need(nrow(running_pubmed_dup()) != 0, '0 articles associated with duplications found.')
      )
      
      tmp_df <- running_pubmed_dup()
      
    }

    if(input$only_omim == 'Yes') {

      tmp_df <- tmp_df %>% 
        select(pmid, everything()) %>%
        left_join(omim_assoc(), by = c('pmid' = 'pubmed_id')) %>%
        filter(omim_assoc != '') %>%
        mutate(omim_assoc = paste0("<a href='", paste0('https://www.omim.org/entry/', omim_assoc),"' target='_blank'>", omim_assoc,"</a>")) %>%
        mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>"))

      vector_colnames <- c('PMID', 'Title','First author', 'Last author', 'N°cites','Journal', 'Published date', 'OMIM entries')
      
    } else {
      
      tmp_df <- tmp_df %>% 
        select(pmid, everything()) %>%
        mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>"))
      
      
      vector_colnames <- c('PMID', 'Title','First author', 'Last author', 'N°cites','Journal', 'Published date')
      
    }
    
    
    datatable(tmp_df, 
              extension = 'Scroller',
              rownames = FALSE, 
              filter = 'top', 
              selection = 'single', 
              escape = FALSE,
              colnames = vector_colnames,
              options = list(
                # autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE,
                scrollY = 200,
                pageLength = 20, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                selection = 'single'
                # columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
    
  })
  
  
  
  
  running_pubtator <- reactive({
    
    req(running_pubmed_del())
    req(running_pubmed_del())
    
    if (input$select_del_dup_frequency == 'deletions') {
      
      input_pubmed <- running_pubmed_del()
        
    } else {
      
      input_pubmed <- running_pubmed_dup()
        
    }
    
    n_pubmed <- nrow(input_pubmed)
    limit_n_articles <- ifelse(n_pubmed > 100, 100, n_pubmed)
    
    tryCatch(
      pubtator_list <- find_pubtator(pmid = input_pubmed$pmid[1:limit_n_articles], bioconcept = 'all', raw_abstract = TRUE),
      error= function(e) stop("The server is not available, please try later.")
    )

    return(pubtator_list)
  })
  
  
  
  output$pubtator_plot_disease <- renderPlot({
    
    entities_df <- running_pubtator() %>% map_dfr(~.[['dataframe']])
    
    p1 <- entities_df %>%
      filter(category == 'Gene') %>%
      distinct(id, word) %>%
      count(word) %>%
      arrange(desc(n)) %>%
      slice(1:10) %>%
      ggplot(aes(reorder(word, n), n)) +
      geom_col(aes( fill = n),color = 'black', show.legend = FALSE) +
      scale_fill_viridis_c() +
      coord_flip() +
      theme_minimal() +
      labs(title = 'Frequency of genetic entities in titles and abstracts',
           x = 'Genetic entity', y = 'Number of articles')
    
    
    p2 <- entities_df %>%
      filter(category == 'Disease') %>%
      distinct(id, word) %>%
      count(word) %>%
      arrange(desc(n)) %>%
      slice(1:10) %>%
      ggplot(aes(reorder(word, n), n)) +
      geom_col(aes( fill = n),color = 'black', show.legend = FALSE) +
      scale_fill_viridis_c() +
      coord_flip() +
      theme_minimal() +
      labs(title = 'Frequency of phenotypic entities in titles and abstracts',
           x = 'Disease entity', y = 'Number of articles')
    
    df_gene <- entities_df %>%
      filter(category == 'Gene') %>%
      distinct(id, word) %>%
      rename(genetic = word) 
    
    df_disease <- entities_df %>%
      filter(category == 'Disease') %>%
      distinct(id, word) %>%
      rename(phenotypic = word)
    
    p3 <- df_gene %>%
      left_join(df_disease, by = 'id') %>%
      count(genetic, phenotypic, sort = TRUE) %>%
      rowwise() %>%
      mutate(tag = paste(genetic, phenotypic, sep = ' - ')) %>%
      ungroup() %>%
      slice(1:10) %>%
      ggplot(aes(reorder(tag, n), n)) +
      geom_col(aes( fill = n), color = 'black', show.legend = FALSE) +
      scale_fill_viridis_c() +
      coord_flip() +
      theme_minimal() +
      labs(title = 'Associations between genetic and phenotypic entities',
           x = 'Entity association', y = 'Number of articles') +
      theme(plot.title = element_text(size=10))
    
    p1 + p2 + p3

  })


  output$plot_net_pubmed <- renderPlot({

    req(isTRUE(input$enable_net))

    if (input$choose_net_source == 'both') {

      tmp_del <- running_pubmed_del() %>% mutate(type = 'deletion')
      tmp_dup <- running_pubmed_dup() %>% mutate(type = 'duplication')

      tmp_to_plot <- bind_rows(tmp_del, tmp_dup)

    } else if (input$choose_net_source == 'deletion') {

      tmp_to_plot <- running_pubmed_del() %>% mutate(type = 'deletion')


    } else {

      validate(
        need(nrow(running_pubmed_dup()) != 0, '0 articles associated with duplications found.')
      )

      tmp_to_plot <- running_pubmed_dup() %>% mutate(type = 'duplication')

    }
    
    tmp_result <- tmp_to_plot %>%
      select(pmid, title) %>%
      unnest_tokens(word, title)
    
    if (input$entity_title_abstract != 'Title') {

      tmp_result <- running_pubtator() %>%
        map_dfr(~.[['dataframe']]) %>%
        rename(pmid = id) %>%
        select(pmid, word) %>%
        bind_rows(tmp_result)

    }
    
    tmp_result <- tmp_result %>%
      anti_join(stop_words, by = 'word') %>%
      pairwise_count(word, pmid, sort = TRUE) %>%
      filter(!str_detect(item1, 'patient'), !str_detect(item2, 'patient')) %>%
      filter(n >= input$min_threshold_cooccurrence)

    validate(
      need(nrow(tmp_result) != 0, 'Please, reduce the minimum threshold of co-occurrence.')
    )
    
    tmp_result %>%
      graph_from_data_frame() %>%
      ggraph(layout = "fr") +
      geom_edge_link(aes(edge_alpha = n, edge_width = n), edge_colour = "cyan4") +
      geom_node_point(size = 5) +
      geom_node_text(aes(label = name), repel = TRUE,
                     point.padding = unit(0.2, "lines")) +
      theme_void()

  })

  abstract_pubmed <- reactive({
    
      validate(
        need(input$del_dup_pubmed_rows_selected != '', "Please, select an article from the table.")
      )
    
    if (input$select_del_dup == 'deletions') {
      
      paper_selected <- running_pubmed_del() %>% 
        slice(input$del_dup_pubmed_rows_selected) %>%
        pull(pmid)
      
    } else {
      
      paper_selected <- running_pubmed_dup() %>% slice(input$del_dup_pubmed_rows_selected) %>%
        pull(pmid)
      
    }
    

    result <- find_pubtator(pmid = paper_selected, bioconcept = 'all')
    
    
    
  })
  
  
  
  output$abstract_html <- renderUI({
    
    HTML(abstract_pubmed()[[1]]$abstract_tagged)

  })
  
  output$legend_html <- renderUI({
    
    # .domain(["CNV", "Enhancer", "lncRNAs", "miRNAs", "TFs"])
    # .range(["#66C2A5","#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"]);'

HTML('<center>
     <span style="color:#66C2A5"> <b> CNV </b> </span> <br>
     <span style="color:#FC8D62"><b>Enhancer</b></span> <br>
     <span style="color:#8DA0CB"><b>lncRNAs</b></span> <br>
     <span style="color:#E78AC3"><b>miRNAs</b></span> <br>
     <span style="color:#A6D854"><b>TFs</b></span>  <br>
     <span style="color:#ae7373"><b>TADs</b></span>
     </center>')

  })
  
  output$abstract_df <- renderDT({
    

    df <- abstract_pubmed()[[1]]$dataframe %>% filter(element == 'abstract') %>%
      select(word, category, identifier) %>%
      rename(entity = word) %>%
      distinct()
    
    datatable(df, 
              colnames = c('Entity', 'Category', 'Identifier'),
              rownames = FALSE, options = list(dom = 't'), class = 'cell-border stripe')
    
  })
    
  
  # output$abstract_pubmed <- renderDT({
  # 
  #   validate(
  #     need(input$del_dup_pubmed_rows_selected != '', "Please, select an article from the table.")
  #   )
  # 
  # 
  #   if (input$select_del_dup == 'deletions') {
  # 
  #     paper_selected <- running_pubmed_del() %>% slice(input$del_dup_pubmed_rows_selected)
  # 
  #   } else {
  # 
  #     paper_selected <- running_pubmed_dup() %>% slice(input$del_dup_pubmed_rows_selected)
  # 
  #   }
  # 
  # 
  #   fetch.pubmed <- entrez_fetch(db = "pubmed", id = paper_selected %>% pull(pmid), rettype = "xml", parsed = T)
  #   # Extract the Abstracts for the respective IDS.
  #   abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
  #     xmlValue(xmlChildren(x)$Abstract))
  # 
  #   tmp_df <- tibble(Title = paper_selected$title, Abstract = abstracts[[1]]) %>% gather('Category', 'Info')
  #   tmp888 <<- tmp_df
  # 
  #   datatable(tmp_df, rownames = FALSE, colnames = '',
  #             options = list(dom = 't'))
  # 
  # })
  
  
  output$dup_pubmed <- renderDT({
    
    
    tmp_df <- running_pubmed_dup()
    tmp_df <- tmp_df %>% mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>")) 

    datatable(tmp_df, 
              rownames = FALSE, 
              filter = 'top', 
              selection = 'single', 
              escape = FALSE,
              extensions = 'Scroller',
              colnames = c('Title','First author', 'Last author', 'N°cites','Journal', 'Published date', 'PMID' ),
              options = list(
                autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE,
                
                pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                selection = 'single'
                # columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
    
  })

  
  data_selected_prev <- reactive({
    
    
    req(coord_user())
    

        data_raw <- hgcn_genes
    
      
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      data_raw <- data_raw  %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap)
      
      
    } else {
      
      
      data_raw <- data_raw  %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap)
    }

    
    data_raw <- data_raw %>% 
      select(-vg, -ensembl_gene_id)
    
    data_raw <- get_perc_overlap(data_raw, coord_user(),
                                 is_a_gene = TRUE)

    
    return(data_raw)
    
  })
  
  data_selected_enhancers <- reactive({
    
    if (!is.null(input$enhancers_on_off)) {
      if (input$enhancers_on_off) {
        
        # # Genes  mapped in CNV
        # genes_cnv <- data_selected_prev()
        # genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # # Genes NOT mapped in CNV
        # if (is.null(input$df_enhancer_rows_all)) {
        #   enhancers_df <- prev_enhancer()
        # } else {
        #   # enhancers_df <- prev_enhancer()[input$df_enhancer_rows_all,]
        #   enhancers_df <- prev_enhancer()
        #   #
        # }
        # genes_no_cnv <- enhancers_df %>% select(gene) %>% distinct() %>% pull()
        # genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        
        genes_no_cnv <- prev_enhancer() %>% filter(is_mapping == 'No') %>% pull(gene)
        
        # ADAPT IT WHEN ADDING OMIM OR OTHERS!!!
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'Enhancer')
        
        
      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }
    
    # Output
    table_output
  })

  data_selected_mirnas <- reactive({
    
    
    
    if (!is.null(input$mirnas_on_off)) {
      if (input$mirnas_on_off) {
        
        # # Genes  mapped in CNV
        # genes_cnv <- data_selected_prev()
        # genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # # Genes NOT mapped in CNV
        # 
        # if (is.null(input$df_mirna_rows_all)) {
        #   mirnas_df <- mirna_raw()
        # } else {
        #   mirnas_df <- mirna_raw()
        #   # mirnas_df <- mirna_raw()[input$df_mirna_rows_all,]
        # }
        # 
        # 
        # genes_no_cnv <- mirnas_df %>% select(gene_symbol) %>% distinct() %>% pull()
        # genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        
        genes_no_cnv <- mirna_raw() %>% filter(is_mapping == 'No') %>% pull(gene_symbol)
        
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'miRNAs')
        
        
      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }
    
    # Output
    table_output
  })
  
  data_selected_lncrnas <- reactive({
    
    
    
    if (!is.null(input$lncrnas_on_off)) {
      if (input$lncrnas_on_off) {
        
        # # Genes  mapped in CNV
        # genes_cnv <- data_selected_prev()
        # genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # # Genes NOT mapped in CNV
        # 
        # if (is.null(input$lncrna_df_rows_all)) {
        #   lncrnas_df <- lncrna_raw()
        # } else {
        #   lncrnas_df <- lncrna_raw()
        #   # mirnas_df <- mirna_raw()[input$df_mirna_rows_all,]
        # }
        # 
        # 
        # genes_no_cnv <- lncrnas_df %>% select(target_symbol) %>% distinct() %>% pull()
        # genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        
        genes_no_cnv <- lncrna_raw() %>% filter(is_mapping == 'No') %>% pull(target_symbol)
        
        
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'lncRNAs')
        
        
      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }

    table_output
  })
  
  data_selected_tfs <- reactive({
    
    if (!is.null(input$tfs_on_off)) {
      if (input$tfs_on_off) {

        genes_no_cnv <- tf_raw() %>% filter(is_mapping == 'No') %>% pull(target)
        
        
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'TFs')
        
        
      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }
    
    # Output
    table_output
  })
  
  data_selected_tads <- reactive({
    

    if (!is.null(input$tads_on_off)) {
      if (input$tads_on_off) {

        genes_no_cnv <- tads_reactive() %>% 
          filter(boundaries_affected == '1/2') %>% 
          filter(no_mapping != ' - ') %>%
          pull(no_mapping)
        
        if (length(genes_no_cnv) == 2) genes_no_cnv <- paste(genes_no_cnv[1], genes_no_cnv[2])
        
        genes_no_cnv <- (genes_no_cnv %>% str_split(pattern = ', '))[[1]] %>% unique()
       
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'TADs')
        


      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }

    # Output
    table_output
  })
  
  
  
  data_selected <- reactive({

    
    table_output <- data_selected_prev() %>% 
      mutate(source = 'CNV') %>% 
      bind_rows(data_selected_enhancers()) %>%
      bind_rows(data_selected_tads()) %>%
      bind_rows(data_selected_tfs()) %>%
      bind_rows(data_selected_mirnas()) %>%
      bind_rows(data_selected_lncrnas())

    table_output
  })

  running_cnv_syndromes <- reactive({
    

    tmp_df <- syndromes_total %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      replace_na(replace = list(variant_class = '-', phenotypes = '-')) %>%
      select(isca_id, chrom, start, end, syndrome_name, variant_class, phenotypes, source, haplo_description,
             triplo_description) %>%
      get_perc_overlap(coord_user(), is_patho = TRUE)
    
  })
  
  
  output$ui_select_cnv_syndrome <- renderUI({
    
    vector_n_dbs <- running_cnv_syndromes() %>%
      count(source) %>% 
      mutate(new_one = paste0(source, ' (', n, ")")) %>%
      pull(new_one) %>%
      str_replace('decipher', 'DECIPHER') %>% 
      str_replace('clingen', 'ClinGen')

    
    if (length(vector_n_dbs) == 0) {
      
      vector_n_dbs <- c('ClinGen (0)','DECIPHER (0)')
      
    } else if (length(vector_n_dbs) == 1 & str_detect(vector_n_dbs, 'ClinGen')) {
      
      vector_n_dbs <- split(c('clingen'), vector_n_dbs)
      
    } else if (length(vector_n_dbs) == 1 & str_detect(vector_n_dbs, 'DECIPHER')) {
      
      vector_n_dbs <- split(c('decipher'), vector_n_dbs)
      
    } else {
      
      vector_n_dbs <- split(c('clingen', 'decipher'), vector_n_dbs)
      
      
    }

    prettyRadioButtons(
      inputId = "select_cnv_syndrome",
      label = '', 
      choices =  vector_n_dbs,
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
    
    
  })
  
  output$cnv_syndromes <- renderDT({
    
    req(running_cnv_syndromes())
    req(input$select_cnv_syndrome)
    
    validate(
      need(nrow(running_cnv_syndromes()) != 0, "No CNV syndromes found."),
      need(!is.null(running_cnv_syndromes()), "No CNV syndromes found.")
      
    )
    
    
    if (input$select_cnv_syndrome == 'decipher') {
      
      tmp_df <- running_cnv_syndromes() %>% 
        filter(source == 'decipher') %>% 
        select(p_overlap, chrom, start, end, syndrome_name, variant_class, phenotypes)
      
      
      datatable(tmp_df, rownames = FALSE,
                colnames = c('Overlap (%)','Chrom', 'Start', 'End', 'CNV syndrome name', 'Variant class', 'Phenotypes'))
      
    } else {
      
      
      tmp_df <- running_cnv_syndromes() %>% 
        filter(source == 'clingen') %>% 
        select(p_overlap, isca_id, chrom, start, end, syndrome_name, haplo_description, triplo_description) %>%
        mutate(isca_id = paste0("<a href='", paste0('https://search.clinicalgenome.org/kb/gene-dosage/region/', isca_id),"' target='_blank'>", isca_id,"</a>"))

      datatable(tmp_df, rownames = FALSE, escape = FALSE,
                colnames = c('Overlap (%)','ISCA ID', 'Chrom', 'Start', 'End', 'CNV syndrome name', 
                             'Haploinsufficiency evidence', 'Triplosensitivity evidence'))
    }
    
    
    
  })
  
  # observeEvent(input$start_analysis,{
  #   removeUI(selector = "#dgenes")
  # })
  
  running_upset_disease <- reactive({
    

    uspset_df <- data_selected() %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      filter(disease == 'Yes') %>%
      select(-source, -disease) %>%
      select(gene, orphanet, dev, genomics_england, omim, clingen) %>%
      pivot_longer(-gene) %>%
      filter(value == 'Yes') %>%
      rename(term = name) %>%
      mutate(term = str_replace(term, 'orphanet', 'ORPHANET'),
             term = str_replace(term, 'dev', 'DECIPHER'),
             term = str_replace(term, 'genomics_england', 'GENOMICS ENGLAND'),
             term = str_replace(term, 'omim', 'OMIM'),
             term = str_replace(term, 'clingen', 'CLINGEN'))
    
    uspset_df
  })
  
  
  output$plot_upset_disease <- renderPlot({
    

    get_upset(running_upset_disease())
    
    
    
  })
  
  
  running_upset_disease_reg <- reactive({
    
    
    uspset_df <- data_selected() %>% 
      select(-start, -end, -chrom) %>%
      filter(source != 'CNV') %>%
      filter(disease == 'Yes') %>%
      select(-source, -disease) %>%
      select(gene, orphanet, dev, genomics_england, omim, clingen) %>%
      pivot_longer(-gene) %>%
      filter(value == 'Yes') %>%
      rename(term = name) %>%
      mutate(term = str_replace(term, 'orphanet', 'ORPHANET'),
             term = str_replace(term, 'dev', 'DECIPHER'),
             term = str_replace(term, 'genomics_england', 'GENOMICS ENGLAND'),
             term = str_replace(term, 'omim', 'OMIM'),
             term = str_replace(term, 'clingen', 'CLINGEN'))
    
    uspset_df
  })
  
  output$plot_upset_disease_reg <- renderPlot({
    

    get_upset(running_upset_disease_reg())
    
    
    
  })
  
  
  
  output$dgenes <- renderDT({
    
    server <- TRUE

    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      select(-source) %>%
      select(band, gene, disease, orphanet, dev, clingen, omim, gwas, p_overlap) %>%
      filter(disease == 'Yes')
    
    data_tmp <- data_input %>% 
      left_join(running_upset_disease() %>% 
                  count(gene, value), by = 'gene') %>%
      select(gene, n, p_overlap)
    
    
    validate(
      need(nrow(data_tmp) > 0, "0 disease genes.")
    )
    
    tmp_output <- datatable(data_tmp, rownames = FALSE, colnames = c('Gene', 'Nº evidences', 'Overlap (%)'),
                            filter = list(position = 'top'), 
                            selection = 'single',
                            options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    
    tmp_output

  })
  
  output$dgenes_reg <- renderDT({
    
    server <- TRUE
    
    
    
    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source != 'CNV') %>%
      select(source, band, gene, disease, orphanet, dev, clingen, omim, gwas) %>%
      filter(disease == 'Yes')

    
    data_tmp <- data_input %>% 
      left_join(running_upset_disease_reg() %>% 
                  count(gene, value), by = 'gene') %>%
      select(source, gene, n)
    
    
    validate(
      need(nrow(data_tmp) > 0, "0 disease genes.")
    )
    
    tmp_output <- datatable(data_tmp, rownames = FALSE, colnames = c('Source','Gene', 'Nº evidences'),
                            filter = list(position = 'top'), 
                            selection = 'single',
                            options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))
    
    tmp_output
    
  })
  
  
  
  output$input_source <- renderUI({
    
    validate(
      need(input$dgenes_rows_selected != '', '')
    )
    
    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      select(-source) %>%
      filter(disease == 'Yes') %>%
      select(gene, dev, omim, orphanet, genomics_england, clingen) %>%
      slice(input$dgenes_rows_selected)
    

    options_sources <- data_input %>% pivot_longer(cols = -gene) %>% filter(value == 'Yes') %>%
      pull(name)
    
    
    give_options <- split(options_sources, 
                  toupper(options_sources) %>% str_replace('_', ' ') %>% str_replace('DEV', 'DECIPHER') %>%
                    str_replace('GENOMICS ENGLAND', 'GENOMICS ENGLAND PanelApp'))
    
    pickerInput(
      inputId = "select_source",
      label = "", 
      choices = give_options
    )
  })
  
  output$input_source_reg <- renderUI({
    
    validate(
      need(input$dgenes_reg_rows_selected != '', '')
    )
    
    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source != 'CNV') %>%
      select(-source) %>%
      filter(disease == 'Yes') %>%
      select(gene, dev, omim, orphanet, genomics_england, clingen) %>%
      slice(input$dgenes_reg_rows_selected)
    
    
    options_sources <- data_input %>% pivot_longer(cols = -gene) %>% filter(value == 'Yes') %>%
      pull(name)
    
    
    give_options <- split(options_sources, 
                  toupper(options_sources) %>% str_replace('_', ' ') %>% str_replace('DEV', 'DECIPHER') %>%
                    str_replace('GENOMICS ENGLAND', 'GENOMICS ENGLAND PanelApp'))
    
    pickerInput(
      inputId = "select_source_reg",
      label = "", 
      choices = give_options
    )
  })
  
  
  output$select_gene_disease <- renderDT({
    
    
    validate(
      need(input$dgenes_rows_selected != '', 'Please, select a disease gene on the left panel.'),
      need(!is.null(input$dgenes_rows_selected), 'No disease target-genes found.'),
      need(!is.null(input$select_source), 'Please, select a evidence source.')
      
    )
    
    req(input$dgenes_rows_selected)
    
    
    gene_selected <- data_selected() %>% 
      filter(source == 'CNV') %>%
      filter(disease == 'Yes') %>%
      slice(input$dgenes_rows_selected) %>%
      pull(gene)
    
    
    if (input$select_source == 'omim') {
      
      tmp_df <- omim %>% 
        replace_na(list(gene_inheritance_mode = '-')) %>%
        filter(gene == gene_selected) %>%
        select(MIM_gene_number, gene_inheritance_mode, MIM_pheno_number, Phenotype) %>%
        mutate(Phenotype = str_remove(Phenotype, '[0-9]{6}')) %>%
        mutate(MIM_gene_number =  paste0("<a href='", 
                                         paste0('https://www.omim.org/entry/', MIM_gene_number),"' target='_blank'>", MIM_gene_number,"</a>")) %>%
        mutate(MIM_pheno_number =  paste0("<a href='", 
                                          paste0('https://www.omim.org/entry/', MIM_pheno_number),"' target='_blank'>", MIM_pheno_number,"</a>")) 
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene MIM number', 'Inheritance', 'Disease MIM number', 'Phenotype'))
      
    } else if (input$select_source == 'orphanet') {
      
      tmp_df <- orphanet_raw %>% 
        filter(gene  == gene_selected) %>%
        replace_na(list(SourceOfValidation = '-')) %>%
        mutate(SourceOfValidation = str_replace_all(SourceOfValidation, '_', '<br>')) %>%
        select(gene, OrphaNumber5, Name6, Name, OrphaNumber, SourceOfValidation ) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=', OrphaNumber5),"' target='_blank'>", gene,"</a>")) %>%
        mutate(OrphaNumber = paste0("<a href='", 
                                    paste0('https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=', OrphaNumber),"' target='_blank'>", OrphaNumber,"</a>") )
      
      
      
      datatable(tmp_df %>% select(-OrphaNumber5), escape = FALSE, rownames = FALSE, colnames = c('Gene',
                                                                       'Description', 
                                                                       # 'ID Orpha gene', 
                                                                       'Disease', 
                                                                       'ID Orpha disease', 
                                                                       'Source of validation'))
      
    }  else if (input$select_source == 'dev') {
      
      tmp_df <- dev_raw %>%
        filter(gene == gene_selected) %>%
        mutate(`gene mim` =  paste0("<a href='", 
                                    paste0('https://www.omim.org/entry/', `gene mim`),"' target='_blank'>", `gene mim`,"</a>")) %>%
        mutate(`organ specificity list` = str_replace_all(`organ specificity list`, ';', '<br>')) %>%
        mutate(pmids = str_replace_all(pmids, ';', '<br>'))
      
      
      datatable(tmp_df, rownames = FALSE, escape = FALSE, colnames = c('Gene', 'Gene MIM number', 'Disease', 'Allelic requirement', 
                                                                       'Organ specificity', 'PMIDS'))
      
    } else if (input$select_source == 'genomics_england') {
      
      tmp_df <- panel_total %>%
        filter(gene == gene_selected) %>%
        mutate(source = str_remove(source, 'gene_panel_')) %>%
        mutate(Phenotypes = str_replace_all(Phenotypes, ';', '<br>')) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://panelapp.genomicsengland.co.uk/panels/entities/', gene),"' target='_blank'>", gene,"</a>")) %>%
        mutate(source = paste0("<a href='", 
                               paste0('https://panelapp.genomicsengland.co.uk/panels/', source),"' target='_blank'>", source,"</a>")) %>%
        mutate(Level4 = paste0(Level4, ' (', source, ')')) %>% 
        select(-source)
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'Gene panel', 'Phenotype associated'))
      
    } else if (input$select_source == 'clingen') {
      
      tmp_df <- data_selected() %>%
        filter(gene == gene_selected) %>%
        select(gene, haplo, triplo) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym=', gene),"' target='_blank'>", gene,"</a>"))
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'Haploinsufficient', 'Triplosensitivity'))
    }
    
  })
  
  output$select_gene_disease_reg <- renderDT({
    
    
    validate(
      need(input$dgenes_reg_rows_selected != '', 'Please, select a disease gene on the left panel.'),
      need(!is.null(input$dgenes_reg_rows_selected), '0 disease genes found.'),
      need(!is.null(input$select_source_reg), '0 disease genes found.')
      
    )
    
    req(input$dgenes_reg_rows_selected)

    gene_selected <- data_selected() %>% 
      filter(source != 'CNV') %>%
      filter(disease == 'Yes') %>%
      slice(input$dgenes_reg_rows_selected) %>%
      pull(gene)
    
    
    if (input$select_source_reg == 'omim') {
      
      tmp_df <- omim %>% 
        replace_na(list(gene_inheritance_mode = '-')) %>%
        filter(gene == gene_selected) %>%
        select(MIM_gene_number, gene_inheritance_mode, MIM_pheno_number, Phenotype) %>%
        mutate(Phenotype = str_remove(Phenotype, '[0-9]{6}')) %>%
        mutate(MIM_gene_number =  paste0("<a href='", 
                                         paste0('https://www.omim.org/entry/', MIM_gene_number),"' target='_blank'>", MIM_gene_number,"</a>")) %>%
        mutate(MIM_pheno_number =  paste0("<a href='", 
                                          paste0('https://www.omim.org/entry/', MIM_pheno_number),"' target='_blank'>", MIM_pheno_number,"</a>")) 
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene MIM number', 'Inheritance', 'Disease MIM number', 'Phenotype'))
      
    } else if (input$select_source_reg == 'orphanet') {
      
      tmp_df <- orphanet_raw %>% 
        filter(gene  == gene_selected) %>%
        replace_na(list(SourceOfValidation = '-')) %>%
        mutate(SourceOfValidation = str_replace_all(SourceOfValidation, '_', '<br>')) %>%
        select(gene, OrphaNumber5, Name6, Name, OrphaNumber, SourceOfValidation ) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=', OrphaNumber5),"' target='_blank'>", gene,"</a>")) %>%
        mutate(OrphaNumber = paste0("<a href='", 
                                    paste0('https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=', OrphaNumber),"' target='_blank'>", OrphaNumber,"</a>") )
      
      
      
      datatable(tmp_df %>% select(-OrphaNumber5), escape = FALSE, rownames = FALSE, colnames = c('Gene',
                                                                       'Description', 
                                                                       # 'ID Orpha gene', 
                                                                       'Disease', 
                                                                       'ID Orpha disease', 
                                                                       'Source of validation'))
      
    }  else if (input$select_source_reg == 'dev') {
      
      tmp_df <- dev_raw %>%
        filter(gene == gene_selected) %>%
        mutate(`gene mim` =  paste0("<a href='", 
                                    paste0('https://www.omim.org/entry/', `gene mim`),"' target='_blank'>", `gene mim`,"</a>")) %>%
        mutate(`organ specificity list` = str_replace_all(`organ specificity list`, ';', '<br>')) %>%
        mutate(pmids = str_replace_all(pmids, ';', '<br>'))
      
      
      datatable(tmp_df, rownames = FALSE, escape = FALSE, colnames = c('Gene', 'Gene MIM number', 'Disease', 'Allelic requirement', 
                                                                       'Organ specificity', 'PMIDS'))
      
    } else if (input$select_source_reg == 'genomics_england') {
      
      tmp_df <- panel_total %>%
        filter(gene == gene_selected) %>%
        mutate(source = str_remove(source, 'gene_panel_')) %>%
        mutate(Phenotypes = str_replace_all(Phenotypes, ';', '<br>')) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://panelapp.genomicsengland.co.uk/panels/entities/', gene),"' target='_blank'>", gene,"</a>")) %>%
        mutate(source = paste0("<a href='", 
                               paste0('https://panelapp.genomicsengland.co.uk/panels/', source),"' target='_blank'>", source,"</a>")) %>%
        mutate(Level4 = paste0(Level4, ' (', source, ')')) %>% 
        select(-source)
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'Gene panel', 'Phenotype associated'))
      
    } else if (input$select_source_reg == 'clingen') {
      
      tmp_df <- data_selected() %>%
        filter(gene == gene_selected) %>%
        select(gene, haplo, triplo) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym=', gene),"' target='_blank'>", gene,"</a>"))
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'Haploinsufficient', 'Triplosensitivity'))
      
      
    }
    
  })
  
  output$dgenes_no_disease <- renderDT({
    
    server <- TRUE
    
    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      select(p_overlap, band, gene, disease, fusil, ohnolog, imprinted, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
             ncgerp)
    
    validate(
      need(nrow(data_input) > 0, "0 no disease genes.")
    )
    
    
    tmp_output <- datatable(data_input, 
                            extensions = 'Scroller',
                            options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                            rownames = FALSE, 
                            colnames = c('Overlap(%)', 'Band', 'Gene', 'Disease', 'Essentiality',
                                         'Ohnolog', 'Imprinted',
                                         'pLI', 'RVIS', 'CCR', 'HI', 'GDI', 'SnIPRE', 'ncRVIS',
                                         'ncGERP'),
                            filter = list(position = 'top'), 
                            selection = 'single') %>%
      formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp'), color = styleInterval(94, c('weight', '#ff7f7f'))) %>%
      formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f')))

    
    tmp_output

  })

  
  output$choose_reg_region <- renderUI({
    
    
    tmp_data <- data_selected() %>% 
      filter(source != 'CNV') %>%
      # filter(disease != 'Yes') %>%
      count(source) %>%
      mutate(source = as.character(source)) %>%
      na.omit()
    
    
    if (nrow(tmp_data) == 0) {
      
      NULL
      
      
    } else {
      
      vector_n_dbs <- split(tmp_data$source, paste(tmp_data$source, paste0('(', tmp_data$n, ')')))
      
      prettyRadioButtons(
        inputId = "select_reg_region",
        label = '', 
        choices =  vector_n_dbs,
        inline = TRUE, 
        status = "primary",
        fill = TRUE
      )
      
    }
    
  })
  
  
  
  
  
  output$genes_from_reg_regions <- renderDT({
    
    server <- TRUE

    validate(
      need(!is.null(input$select_reg_region), "0 non-disease target genes found."),
      need(length(input$select_reg_region) != 0, "0 non-disease target genes found.")
    )
    
    data_input <- data_selected() %>% 
      select(-start, -end, -chrom) %>%
      filter(source ==  input$select_reg_region) %>%
      select(band, gene, disease, fusil, ohnolog, imprinted, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
             ncgerp)
    
    validate(
      need(nrow(data_input) > 0, "0 no disease genes.")
    )
    
    
    tmp_output <- datatable(data_input, rownames = FALSE, 
                            extensions = 'Scroller',
                            options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                            colnames = c('Band', 'Gene', 'Disease', 'Essentiality',
                                         'Ohnolog', 'Imprinted',
                                         'pLI', 'RVIS', 'CCR', 'HI', 'GDI', 'SnIPRE', 'ncRVIS',
                                         'ncGERP'),
                            filter = list(position = 'top'), 
                            selection = 'single') %>%
      formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp'), color = styleInterval(94, c('weight', '#ff7f7f'))) %>%
      formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f')))
    
    
    tmp_output

    # server <- TRUE
    # data_input <- data_selected() %>% filter(source ==  input$select_reg_region) %>%
    #   select(-start, -end, -chrom) %>%
    #   select(band, gene, disease, essent, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
    #          ncgerp) %>%
    #   filter(disease == 'No')
    # 
    # datatable(data_input, rownames = FALSE, filter = list(position = 'top'),
    #           colnames = c('Band', 'Gene', 'Disease', 'Essential', 'pLI', 'RVIS', 'CCR', 'HI', 'GDI', 'SnIPRE', 'ncRVIS',
    #                        'ncGERP')) %>%
    #   formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp'), color = styleInterval(94, c('weight', '#ff7f7f'))) %>%
    #   formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f'))) %>%
    #   formatStyle(c('disease'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
    
    # formatStyle(c('disease', 'haplo', 'triplo', 'omim', 'dev', 'fda', 'gwas'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
    
    # 
    #   }
    # }
    # )
  })
  
  # output$genes_from_tads <- renderDT({
  #   
  #   if (!is.null(input$tads_on_off)) {
  #     if (input$tads_on_off) {
  #       
  # 
  #       server <- TRUE
  #       data_input <- data_selected() %>% filter(source == 'TAD') %>%
  #         select(-start_position, -end_position, -chrom) %>%
  #         select(band, gene, disease, haplo, triplo, dev, fda, omim, fda, gwas, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
  #                ncgerp, p_overlap)
  #       
  #       datatable(data_input, rownames = FALSE, filter = 'top',
  #                 selection = 'single',
  #                 options = list(
  #                   pageLength = 5,
  #                   # autoWidth = TRUE,
  #                   style = 'bootstrap',
  #                   list(searchHighlight = TRUE),
  #                   stateSave = FALSE
  #                   # colnames = c('Entrez id', 'Band', 'Gene', 'pLI', 'Database', 'CNV size', 'Percentage Overlap (%)')
  #                   # columnDefs = list(list(className = 'dt-center', targets = '_all'))
  #                 )) %>%
  #         formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp'), color = styleInterval(95, c('weight', '#ff7f7f'))) %>%
  #         formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f'))) %>%
  #         formatStyle(c('disease', 'haplo', 'triplo', 'omim', 'dev', 'fda', 'gwas'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
  #     }
  #   }
  #   # )
  # })
  
  output$score_references <- renderDT({
    
    datatable(ref_scores, rownames = FALSE,
              options = list(
                pageLength = 5, autoWidth = TRUE, style = 'bootstrap',
                selection = 'single',
                columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
  })
  

  
  output$plot_size <- renderPlot({
    
    req(coord_user())
   
    
    size_cnv_query = coord_user() %>% mutate(length_cnv_input = end - start + 1)

    if (input$select_density == 'global') {

      ridges_home +
        geom_vline(data = size_cnv_query,
                   aes(xintercept = length_cnv_input), linetype = 2, color = 'red', size = 1.5)
      
      
    } else {
      check_cnv_df() %>%
        mutate(source =
                 case_when(
                   source == 'dgv' ~ 'DGV',
                   source == 'decipher' ~ 'DECIPHER',
                   source == 'gnomad_v2.1' ~ 'gnomAD v2.1',
                   source == 'decipher_control' ~ 'DECIPHER Control'
                 )) %>%
        bind_rows(clinvar_variants %>% 
                    filter(length_cnv >= 50) %>% 
                    bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
                    mutate(source = paste('ClinVar','-', clinical)) %>% 
                    select(source, length_cnv)) %>%
        ggplot(aes(length_cnv, y = source)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = source), alpha = 0.6, 
                            show.legend = FALSE, size = 1.25, bandwidth = 0.304) +
        geom_vline(data = size_cnv_query, aes(xintercept = length_cnv_input), linetype = 2, color = 'red', size = 1.5) +
        scale_x_log10() +
        scale_y_discrete(expand = c(0.01, 0)) +
        scale_fill_viridis_d() +
        xlab('log10(CNVs size)') +
        ylab('Database') +
        theme_ridges()
      
    }
    
  })
  
  output$tissue_gtex <- renderPlot({
    
    req(!is.null(input$gtex_gene_tissue))

    if (input$gtex_gene_tissue == 'Gene' ) {
      
      req(!is.null(input$input_gtex_gene))

      filtered_gene <- input$input_gtex_gene
      
      p <- gtex %>%
        filter(gene == filtered_gene) %>%
        ggplot(aes(reorder(tissue, -value), value)) +
        geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
        theme_fancy() +
        xlab('Tissue') +
        ylab(paste('Median TPM')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position='none') +
        ggtitle(paste('Gene expression: ', filtered_gene))
      
    } else {
      
      req(!is.null(input$input_gtex_tissue))
      
      filtered_genes_cnv <- data_selected() %>% pull(gene)
      filtered_tissue <- input$input_gtex_tissue
      
      p <- gtex %>%
        filter(tissue == filtered_tissue) %>%
        filter(gene %in% filtered_genes_cnv) %>%
        ggplot(aes(reorder(gene, -value), value)) +
        geom_col(aes(fill = tissue), color = 'black', show.legend = FALSE) +
        theme_fancy() +
        xlab('Tissue') +
        ylab(paste('Median TPM')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position='none') +
        ggtitle(paste('Tissue: ', filtered_tissue))
      
    }
    
    p
    
  })
  
  running_tsea <- reactive({
    
    req(isTRUE(input$enable_tsea))
    
    
    vector_genes <- data_selected() %>% pull(gene)
    vector_background <- hgcn_genes %>% pull(gene)
    
    gs<- GeneSet(geneIds= vector_genes,
                organism="Homo Sapiens",
                geneIdType=SymbolIdentifier())
    
    background_genes <- GeneSet(geneIds= vector_background,
                 organism="Homo Sapiens",
                 geneIdType=SymbolIdentifier())
    
    output <- teEnrichment(inputGenes = gs, 
                           rnaSeqDataset = if_else(input$tissue_expression_dbs == 'gtex', 1, 3),
                           backgroundGenes = background_genes)
    
    
    seEnrichmentOutput<-output[[1]]
    
    enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
    
    enrichmentOutput$Tissue<-row.names(enrichmentOutput)
    
    enrichmentOutput
    
    
  })
  
  output$plot_tsea <- renderPlot({
    
    
    validate(
      need(nrow(running_tsea() %>% filter(Log10PValue > 1.301)) != 0, "0 enriched tissues found.")
    )
    
    
    
    running_tsea()  %>%
      ggplot(aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
      geom_bar(stat = 'identity', color = 'black')+
      labs(x='', y = '-LOG10(P-Adjusted)')+
      theme_fancy()+
      theme(legend.position="none")+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
    
  })
  
  
  running_hpa <- reactive({
    
    vector_genes <- data_selected() %>% pull(gene)
    
    tmp_df <- hpa %>% filter(gene %in% vector_genes)
    
    tmp_df

  })
  
  output$tissue_hpa <- renderDT({
    
    tmp_df <- running_hpa()
    
    if (input$gene_yes_no == 'Yes') {
      
      req(input$input_gene_tissue)
      filtered_gene <- input$input_gene_tissue
      tmp_df <- tmp_df %>% filter(gene == !!filtered_gene)
    }
    
    if (input$tissue_yes_no == 'Yes') {
      
      req(input$input_tissue)
      filtered_tissue <- input$input_tissue
      tmp_df <- tmp_df %>% filter(tissue == !!filtered_tissue)
    }
    
    validate(
      need(nrow(tmp_df) != 0, "No data found.")
    )


    datatable(tmp_df, 
              extensions = 'Scroller',
              options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
              colnames = c('Gene','Tissue', 'Cell type', 'Level', 'Reliability'),
              filter = 'top')
    
    
    
  })
  

  
  
  
  
  output$model_genes <- renderDT({
    
    validate(
      need(nrow(model_genes_phenotype()) != 0, '0 genes found.')
    )
    
    
    df_tmp <- model_genes_phenotype() %>% 
      select(-entrez_id, -gene_mouse) %>%
      mutate(mgi = paste0("<a href='", paste0('http://www.informatics.jax.org/marker/', mgi),"' target='_blank'>", mgi,"</a>")) %>%
      mutate(term = paste0("<a href='", paste0('http://www.informatics.jax.org/vocab/mp_ontology/', term),"' target='_blank'>", term,"</a>"))
    
    
    datatable(df_tmp, 
              extensions = 'Scroller',
              escape = FALSE, 
              rownames = FALSE, 
              colnames = c('Source','Human ortholog gene', 'Mouse gene', 'Mouse phenotype id', 'Phenotype description'),
              options = list(
                deferRender = TRUE, scrollY = 200, scroller = TRUE,
                columnDefs = list(list(className = 'dt-center', targets = 0:4) )))

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

  output$data <- renderTable({
    mtcars[, c("mpg", input$variable), drop = FALSE]
  }, rownames = TRUE)
  
  
  
  plot_chrom_react <- reactive({
    
  req(input$input_geno_karyo != 'Multiple coordinates (NGS)')
  req(!is.null(coord_user()))
  req(nrow(coord_user()) == 1)

    start_coordinates <- coord_user() %>% pull(start)
    end_coordinates <- coord_user() %>% pull(end)
    chrom_coordinates <- coord_user() %>% pull(chrom)
    chrom_coordinates <- paste0('chr', chrom_coordinates)

    ideoTrack <- Gviz::IdeogramTrack(genome= "hg19", chromosome = chrom_coordinates)
    Gviz::plotTracks(ideoTrack, from= start_coordinates , to = end_coordinates, 
                     showBandId=TRUE, cex.bands = 0.5)
    
  })
  
  
  output$plot_chrom <- renderPlot({
    
    plot_chrom_react()
    
    
  })



  
  
 
  
  output$choose_geno_karyo1 <- renderUI({
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      textInput(
        inputId = "int_start",
        label = "Genomic interval - Start",
        value = '153613593')
      
    } else {
      
      karyotype_filtered <- as.list(coord_cytobands %>% filter(Chrom == input$input_chrom) %>% select(Name))
      
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

    tablerInfoCard(
      width = 12,
      value =  paste(nrow(data_selected()), 'Genes'),
      status = "primary",
      icon = "database",
      description = 'Number of genes - CNV')

  })

  
  output$genes_enhancers_selected <- renderUI({

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
      status = "primary",
      icon = "database",
      # description =  name_region
      description = 'Total number of genes'
    )
    
    
  })
  
  
  output$ref_user_genes_cnv <- renderUI({
    
    
    rows_selected <- input$dgenes_rows_all

    
    tmp_df <- data_selected() %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV')
    
    
    if (length(rows_selected) == nrow(tmp_df) | is.null(rows_selected)) {
      tablerInfoCard(
        width = 12,
        value = paste0(nrow(tmp_df), " genes"),
        status = "success",
        icon = "database",
        description = 'Genes found in CNV(s)'
      )
    } else {
      tablerInfoCard(
        width = 12,
        value = paste0(nrow(tmp_df), " genes"),
        status = "success",
        icon = "database",
        description = 'Genes found in CNV(s)'
      )
    }
    
  })
  
  
  output$n_variants <- renderUI({
    
    tablerInfoCard(
      width = 12,
      value = paste0(nrow(cnv_file_to_analyze()), " variants"),
      status = "success",
      icon = "database",
      description = 'Nº of variants selected'
    )
  })
  
  output$checking_quality_file <- renderUI({
    
    original_n_rows <- reading_cnv_file() %>% nrow()
    
    dup_n_rows <-  original_n_rows - (reading_cnv_file() %>% distinct() %>% nrow())
    
    if (dup_n_rows > 0) {

      tablerInfoCard(
        width = 12,
        value = paste0(dup_n_rows, " duplicated"),
        status = "success",
        icon = "database",
        description = 'Nº of rows duplicated'
      )
      
      }
  })
  
  
  generate_cnvscore <- reactive({
    
    req(input$input_geno_karyo != 'Multiple coordinates (NGS)')
    
    tmp_start <- coord_user() %>% pull(start)
    tmp_end <- coord_user() %>% pull(end)
    tmp_chrom <- coord_user() %>% pull(chrom)
    
    test10 <<- tmp_start
    test11 <<- tmp_end
    test12 <<- tmp_chrom
    
    api_result <- POST('http://3.68.213.5:3838/classifier',
                       body = paste0('input_chrom=', tmp_chrom,
                                     '&input_start=', as.integer(tmp_start), 
                                     '&input_end=',  as.integer(tmp_end), 
                                     '&input_type=deletion'))
    
    api_result <- as_tibble(content(api_result)[[1]])

  })
  
  output$ref_user_cnvscore <- renderUI({
    
    test14 <<- generate_cnvscore()
    
    tmp_df <- generate_cnvscore()
    
    tablerInfoCard(
      width = 12,
      value =  paste0(ifelse(tmp_df$.pred_pathogenic >= 0.5, 'Pathogenic CNV \n', 'Benign CNV \n'),
                      'Uncertainty:', tmp_df$sd),
      status = "primary",
      icon = "database",
      description =  paste(tmp_df$chrom, tmp_df$start, tmp_df$end)
    )

  })
  
  output$cnvscore_rules <- renderDataTable({
    
    tmp_df <- generate_cnvscore()
    
    tmp_df <- str_split(tmp_df$rules, ', ')[[1]] %>% as_tibble()
    
    datatable(tmp_df)

    
  })
  
  output$ref_user_region <- renderUI({
    
    req(input$input_geno_karyo != 'Multiple coordinates (NGS)')
    
    tmp_start <- coord_user() %>% pull(start)
    tmp_end <- coord_user() %>% pull(end)
    tmp_chrom <- coord_user() %>% pull(chrom)
    
    
    tablerInfoCard(
      width = 12,
      value =  paste0(tmp_chrom,':', tmp_start, '-', tmp_end),
      status = "primary",
      icon = "database",
      description =  'Current region displayed'
    )
    
    
    
  })
  
  
  output$ref_user_region_file <- renderUI({
    
    tmp_start <- coord_user() %>% pull(start)
    tmp_end <- coord_user() %>% pull(end)
    tmp_chrom <- coord_user() %>% pull(chrom)
    
    tablerInfoCard(
      width = 12,
      value =  coord_user() %>% nrow(),
      status = "primary",
      icon = "database",
      description =  'Nº region(s) displayed'
    )
  })
  
  
  output$ref_user_length <- renderUI({
    
    
    if (nrow(coord_user()) > 1) {
      
      length_region <- coord_user() %>% 
        mutate(length_region = end - start + 1) %>% pull(length_region) %>% sum()
      
      tmp_description <- 'Total length of the genomic regions'
      
    } else {
      
      length_region <- coord_user() %>% mutate(length_region = end - start + 1) %>% pull(length_region)
      
      tmp_description <- 'Length of the genomic region'
    }
    
     
    
    if (length_region >= 1e6) {
      length_region <- round(length_region / 1e6, 2)
      tmp_out <- paste(length_region, 'Mb')
    } else {
      length_region <- round(length_region / 1e3, 2)
      tmp_out <- paste(length_region, 'Kbs')
    }
    
    tablerInfoCard(
      width = 12,
      value =  tmp_out,
      status = "primary",
      icon = "database",
      description =  tmp_description
    )
    
    
  })
  
  output$ref_user_cytoband <- renderUI({
    
    req(input$input_geno_karyo != 'Multiple coordinates (NGS)')
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      
      tmp_cyto <- coord_cytobands %>%
        rename(chrom = Chrom, start = Start, end = End) %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        select(Name) %>%
        pull() %>%
        paste(collapse = ', ')
      
      
    } else {
      
      
      tmp_cyto <- input$input_karyotype
      
      
    }
    
    tablerInfoCard(
      width = 12,
      value =  tmp_cyto,
      status = "primary",
      icon = "database",
      description =  'Cytoband(s) mapping'
      
    )
    
    
  })

  
  output$counter_header <- renderUI({

    tmp_tbl <- data_selected() %>% 
      filter(source != 'CNV') %>%
      count(source) %>%
      na.omit()
    

    tmp_n_genes_reg <- tmp_tbl %>% nrow()
    
    req(tmp_n_genes_reg > 0)
    

    tmp_desc_2 <- paste(paste(tmp_tbl$n, 'Target-genes', '-', tmp_tbl$source ), collapse =" <br/> ")

    tablerInfoCard(
      width = 12,
      value = paste0('+', tmp_tbl %>% pull(n) %>% sum(), " genes"),
      status = "warning",
      icon = "database",
      description =  HTML(tmp_desc_2)
    )

  })
  
  output$n_genes_tf_added <- renderUI({
    
    req(nrow(data_selected_tfs()) > 0)
    
    tmp_df <- data_selected_tfs()
    
    tablerInfoCard(
      width = 12,
      value = paste0('+', nrow(tmp_df), " genes"),
      status = "warning",
      icon = "database",
      # description =  name_region
      description = 'Target-genes TFs'
    )
    
  })
  
  
  
  
  output$n_genes_tad_added <- renderUI({
    
    req(nrow(data_selected_tads()) > 0)
    
    tmp_df <- data_selected_tads()
    
    tablerInfoCard(
      width = 12,
      value = paste0('+', nrow(tmp_df), " genes"),
      status = "warning",
      icon = "database",
      # description =  name_region
      description = 'Genes disrupted in TADs'
    )
    
  })
  
  
  mirna_raw <- reactive({
    
    # req(input$start_analysis > 0)

    genes_cnv <- data_selected_prev() %>% pull(gene)
    
    
    data_tmp <- mirtarbase %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      filter(gene_symbol %in% hgcn_genes$gene) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      select(-contains('.source')) %>%
      mutate(is_mapping = if_else(gene_symbol %in% genes_cnv, 'Yes', 'No'))

  })
  
  
  
  output$df_mirna <- renderDT({
    
    validate(
      need(nrow(mirna_raw()) > 0, '0 miRNAs found.')
    )

    tmp_df <- mirna_raw() %>%
      mutate(references = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', references),
                                 "' target='_blank'>", references,"</a>")) %>%
      mutate(name = paste0("<a href='", paste0('http://www.mirbase.org/textsearch.shtml?q=', name),
                           "' target='_blank'>", name,"</a>")) %>%
      mutate(id = paste0("<a href='", paste0('http://mirtarbase.cuhk.edu.cn/php/detail.php?mirtid=', id),
                         "' target='_blank'>", id,"</a>")) %>%
      select(id, name, chrom, start, end, gene_symbol, references, is_mapping, experiment)
    

    datatable(tmp_df, rownames = FALSE, escape = FALSE,
              extensions = 'Scroller',
              options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
              filter = 'top',
              colnames = c('ID', 'Name', 'Chrom', 'Start', 'End',
                           'Target-gene', 'Reference', 'Mapping query?', 'Validation experiment'))
    
    
    
  })
  
  tf_raw <- reactive({
    
    # req(input$start_analysis > 0)
    
    genes_cnv <- data_selected_prev() %>% pull(gene)


    data_tmp <- trrust %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      filter(target %in% hgcn_genes$gene) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      select(-target_chrom, -target_start, -target_end) %>%
      mutate(reference = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', reference),"' target='_blank'>", reference,"</a>")) %>%
      select(-contains('.source')) %>%
      mutate(is_mapping = if_else(target %in% genes_cnv, 'Yes', 'No'))
    

  })
  
  
  output$tf_df <- renderDT({
    
    validate(
      need(nrow(tf_raw()) > 0, '0 Transcription factors (TFs) found.')
    )
    
    datatable(tf_raw(), 
              extensions = 'Scroller',
              filter = 'top',
              options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE),
              rownames = FALSE, escape = FALSE,  
              colnames = c('TF', 'Chrom', 'Start', 'End', 'Target-gene', 'Mechanism', 'Reference', 'Mapping query?'))

  })
  
  
  lncrna_raw <- reactive({
    
    # req(input$start_analysis > 0)
    
    genes_cnv <- data_selected_prev() %>% pull(gene)
    

    data_tmp <- lncrna_target %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      filter(target_symbol %in% hgcn_genes$gene) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      mutate(is_mapping = if_else(target_symbol %in% genes_cnv, 'Yes', 'No')) %>%
      select(-.source)
    
    
    
    
  })
  
  
  
  output$lncrna_df <- renderDT({
    
    validate(
      need(nrow(lncrna_raw()) > 0, '0 lncRNAs found.')
    )

    datatable(lncrna_raw(), 
              extension = 'Scroller',
              rownames = FALSE,
              colnames = c('Symbol', 'Chromosome', 'Start', 'End', 'Ensembl ID', 'Target symbol', 
                           'Tissue origin', 'Disease state', 'PMID', 'Mapping query'),
              options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE,
                pageLength = 5, autoWidth = TRUE, list(searchHighlight = TRUE)))
    
    
    
  })
  
  
  output$n_lncrna <- renderUI({
    
    tablerStatCard(
      value =  nrow(lncrna_raw()),
      title = "lncRNAs",
      trend = NULL,
      width = 12
    )
    
  })
  
  output$n_mirna <- renderUI({
    
    tablerStatCard(
      value =  mirna_raw() %>% nrow(),
      title = "miRNAs",
      trend = NULL,
      width = 12
    )
    
  })
  
  output$n_tfs <- renderUI({
    
    tablerStatCard(
      value =  tf_raw()  %>% nrow(),
      title = "TFs",
      trend = NULL,
      width = 12
    )
    
  })

  
  prev_enhancer <- reactive({
    

    genes_cnv <- data_selected_prev() %>% pull(gene)

    data_tmp <- df_enhancers %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      filter(gene %in% hgcn_genes$gene) %>%
      select(id,chrom, start, end, gene, everything()) %>%
      mutate(is_mapping = if_else(gene %in% genes_cnv, 'Yes', 'No')) %>%
      select(-contains('source'), -contains('delete'),-oe, -.overlap, -score_enh, -score)
    
    data_tmp
  })
  
  output$n_enhancer <- renderUI({
    

    data_tmp <- prev_enhancer() %>% select(id) %>% distinct() %>% pull(id)
    
    
    tablerStatCard(
      value =  length(data_tmp),
      title = "Enhancers",
      # trend = -10,
      width = 12
    )
  })

  
  output$p100_enhancer <- renderPlot({
    
    validate(
      need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
    )
    
    score_filtered <- prev_enhancer() %>% 
      slice(input$df_enhancer_rows_selected) %>% 
      select(phast100) %>% 
      pull()

    plot_p100 +
      theme_fancy() +
      geom_vline(xintercept = score_filtered, color = 'red')
  })
  
  output$p46pla_enhancer <- renderPlot({
    
    validate(
      need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
    )
    
    score_filtered <- prev_enhancer() %>% 
      slice(input$df_enhancer_rows_selected) %>% 
      select(phast46pla) %>% 
      pull()
    
    plot_p46pla +
    theme_fancy() +
      geom_vline(xintercept = score_filtered, color = 'red')


  })
  
  output$p46pri_enhancer <- renderPlot({
    
    validate(
      need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
    )
    
    score_filtered <- prev_enhancer() %>% 
      slice(input$df_enhancer_rows_selected) %>% 
      select(phast46pri) %>% 
      pull()
    
    plot_p46pri +
      xlab('Phast46way primate score') +
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
  
  output$df_enhancer <- renderDT({
    
    validate(
      need(nrow(prev_enhancer()) > 0, '0 enhancers found.')
    )
    

    datatable(prev_enhancer(),
              extensions = 'Scroller',
              rownames = FALSE, filter = 'top', selection = 'single',
              colnames = c('ID Enhancer', 'Chrom',
                           'Start',  'End','Target-gene', 'Phast100way', 'Phast46way Placental', 'Phast46way Primates', 'Mapping query?'),
              options = list(
                deferRender = TRUE, scrollY = 200, scroller = TRUE,
                pageLength = 5, 
                list(searchHighlight = TRUE)
              ))
    
    
  })
  
  number_tads <- reactive({

    # req(input$start_analysis > 0)

    n_tads <- tads_reactive()
    
    if (is.null(n_tads)) {
      n_tads <- 0
    } else {
      n_tads <- nrow(n_tads)
    }
    n_tads

  })
  
  output$n_tads <- renderUI({
    
    
    tablerStatCard(
      value =  number_tads(),
      title = "TADs",
      width = 12
    )

  })
  
  
  output$n_target_enh <- renderUI({
    
    
    n_target_genes <- prev_enhancer() %>% select(gene) %>% distinct() %>% nrow()

    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Unique target-genes of enhancers'
      
    )
    
  })
  
  output$n_target_enh_notmap <- renderUI({
    
    
    n_target_genes <- prev_enhancer() %>% 
      filter(is_mapping == 'No') %>% 
      select(gene) %>% 
      distinct() %>% 
      nrow()
    
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Not mapping the CNV(s)'
      
    )

  })
  
  
  output$n_target_mirna <- renderUI({
    
    
    n_target_genes <- mirna_raw() %>% select(gene_symbol) %>% distinct() %>% nrow()
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Unique target-genes of miRNAs'
      
    )
    
  })
  
  output$n_target_mirna_notmap <- renderUI({
    
    
    n_target_genes <- mirna_raw() %>% filter(is_mapping == 'No') %>% select(gene_symbol) %>% distinct() %>% nrow()
    
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Not mapping with the variant(s)'
      
    )
    
    
  })
  
  
  
  
  output$n_target_tf <- renderUI({
    
    
    n_target_genes <- tf_raw() %>% select(target) %>% distinct() %>% nrow()
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Unique target-genes of TFs'
      
    )
    
  })
  
  output$n_target_tf_notmap <- renderUI({
    
    
    n_target_genes <- tf_raw() %>% filter(is_mapping == 'No') %>% 
      select(target) %>% distinct() %>% nrow()
    
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Not mapping with the variant(s)'
      
    )

    
  })
  
  
  output$n_target_lncrna <- renderUI({
    
    
    n_target_genes <- lncrna_raw() %>% select(target_symbol) %>% distinct() %>% nrow()
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Unique target-genes of lncRNAs'
      
    )
    
  })
  
  
  output$n_target_lncrna_notmap <- renderUI({
    
    n_target_genes <- lncrna_raw() %>% 
      filter(is_mapping == 'No') %>% 
      select(target_symbol) %>% 
      distinct() %>% 
      nrow()

    tablerInfoCard(
      width = 12,
      value =  paste(n_target_genes, 'target-genes'),
      status = "primary",
      icon = "database",
      description =  'Not mapping with the variant(s)'
      
    )
    
    
  })
  
  output$switch_tads <- renderUI({

    switchInput(
      inputId = "tads_on_off",
      label = "Add genes to analysis?",
      inline = TRUE,
      width = 'auto',
      value = FALSE
    )
    
  })
  
  
  output$switch_tfs <- renderUI({
    
    switchInput(
      inputId = "tfs_on_off",
      label = "Add target genes to analysis?",
      inline = TRUE,
      width = 'auto',
      value = FALSE
    )
    
  })
  
  
  output$switch_lncrnas <- renderUI({
    
    switchInput(
      inputId = "lncrnas_on_off",
      label = "Add target genes to analysis?",
      inline = TRUE,
      width = 'auto',
      value = FALSE
    )
    
  })
  output$switch_mirnas <- renderUI({
    
    switchInput(
      inputId = "mirnas_on_off",
      label = "Add target genes to analysis?",
      inline = TRUE,
      width = 'auto',
      value = FALSE
    )
    
  })

  tads_reactive <- reactive({

    n_hits <- coord_user() %>%
      bed_intersect(tad %>%
                      pivot_longer(-c(id, chrom), names_to = 'coord', values_to = 'start') %>%
                      mutate(end = start)) %>%
      select(id.y) %>%
      count(id.y)

    if (length(n_hits$id.y) == 0) return(tibble())
    

    tmp_tads <-   tad %>% 
      filter(id %in% n_hits$id.y) %>%
      rowwise() %>%
      mutate(vector_genes = paste(bed_intersect(hgcn_genes, tibble('chrom' = chrom,
                                                             'start' = start,
                                                             'end' = end)) %>% pull(gene.x),
                                  collapse = ', ')) %>%
      mutate(no_mapping = paste(str_split(vector_genes, pattern = ', ')[[1]][!str_split(vector_genes, pattern = ', ')[[1]] %in%  data_selected_prev()$gene], collapse = ', ')) %>%
      mutate(no_mapping = if_else(no_mapping == '', ' - ', no_mapping)) %>%
      left_join(n_hits %>% rename(boundaries_affected = n), by = c('id' = 'id.y')) %>%
      select(-id) %>%
      mutate(boundaries_affected = paste0(boundaries_affected, '/2'))
    
    
      tmp_tads
  })
  
  output$df_tads <- renderDT({
    
    validate(
      need(nrow(tads_reactive()) > 0, 'No TADs disrupted found.')
    )
    
    datatable(tads_reactive(), 
              # extensions = 'Scroller',
              # options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
              rownames= FALSE,
              colnames = c('Chromosome', 'Start', 'End', 'Genes mapping TAD', 'Genes no mapping CNV(s)',
                           'TAD boundaries disrupted')
    )
    
  })

  
  output$n_clinvar <- renderUI({
    
    
    tablerStatCard(
      value =  nrow(running_clinvar_no_cnv()),
      title =  'ClinVar variants',
      # trend = -10,
      width = 12
    )
    
    
    # tablerInfoCard(
    #   value =  paste(n_clinvar_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Genes found in ClinVar",
    #   width = 12
    # )
    
  })
  
  
  
  output$n_syndromes <- renderUI({
    
    
    n_syndromes <- running_cnv_syndromes()
    
    
    tablerStatCard(
      value =  nrow(n_syndromes),
      title =  'CNV syndromes',
      # trend = -10,
      width = 12
    )

  })
  
  output$n_disease <- renderUI({
    
    
    n_disease_yes <- data_selected() %>% filter(disease == "Yes") %>% nrow()
    n_total <- nrow(data_selected())  
    
    tablerStatCard(
      value =  paste(n_disease_yes, n_total, sep = '/'),
      title =  'Disease genes',
      # trend = -10,
      width = 12
    )
  })
  

  
  output$n_gwas <- renderUI({
    
    
    # n_gwas_yes <- data_selected() %>% filter(gwas == 'Yes') %>% nrow()
    # n_total <- nrow(data_selected())
    
    tablerStatCard(
      # value =  paste(n_gwas_yes, n_total, sep = '/'),
      value = nrow(running_gwas()),
      title =  'GWAS Catalog variants',
      # trend = -10,
      width = 12
    )
    # 
    # tablerInfoCard(
    #   value =  paste(n_gwas_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Genes found in GWAS Catalog",
    #   width = 12
    # )
  })
  
  output$hpo_unique <- renderUI({
    
    
    
    hpo_yes <-  data_selected() %>% select(gene) %>% pull()
    hpo_genes_filter <- hpo_genes %>% 
      filter(gene %in% hpo_yes) %>%
      select(hp) %>%
      distinct()
    
    tablerStatCard(
      value =  nrow(hpo_genes_filter),
      title =  'Unique HPO terms found',
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  nrow(hpo_genes_filter),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Unique HPO terms found",
    #   width = 12
    # )
    
  })
  
  
  
  output$hpo_unique_genes <- renderUI({
    
    hpo_yes <-  data_selected() %>% select(gene) %>% pull()
    
    n_genes <- hpo_genes %>% 
      filter(gene %in% hpo_yes) %>% 
      select(gene) %>% 
      distinct() %>%
      nrow()
    
    tablerStatCard(
      value = n_genes,
      title =  'Genes previously associated with at least 1 HPO term',
      # trend = -10,
      width = 12
    )
  })
  

  
  
  hpo_filter <- reactive({
    
    hpo_yes <-  data_selected() %>% select(gene) %>% pull()
    
    hpo_genes_filter <- hpo_genes %>% 
      filter(gene %in% hpo_yes)
    
    validate(
      need(nrow(hpo_genes_filter) > 0, 'No gene-disease associations found.')
    )
    
    hpo_genes_filter
  })
  
  
  
  
  running_sim_score <- reactive({
    
    req(input$input_inheritance != '')
    
    
      validate(
        need(nrow(data_selected()) != 0, "0 genes found.")

      )

    tmp_tbl <- data_selected() %>% 
      select(source, gene) %>%
      left_join(hpo_genes %>% select(identifier, gene, hp), by = 'gene') %>%
      rename(hp_gene = hp) %>%
      na.omit()
    
    validate(
      need(nrow(tmp_tbl) > 0, 'No gene-disease associations found.')
    )
    
    if (input$input_inheritance != 'Any') {
      
    keep_identifiers <- tmp_tbl %>% select(identifier) %>%
      left_join(hpo_omim %>% select(identifier, hp)) %>%
      filter(hp %in% input$input_inheritance) %>%
      pull(identifier)
    
    tmp_tbl <- tmp_tbl %>% filter(identifier %in% keep_identifiers)
    
    validate(
      need(nrow(tmp_tbl) > 0, 'No gene-disease associations found.')
    )
    
    }

    if (length(input$chosen_hp) == 0) {
      
      return_tbl <- tmp_tbl %>% mutate(sim_gene = NA,
                         sim_mim = NA)
      
    } else {
      
      to_sim_gene <- tmp_tbl %>% 
        select(gene, hp_gene)
      
      to_sim_disease <- tmp_tbl %>% 
        select(identifier) %>%
        left_join(hpo_omim %>% select(identifier, hp)  %>% rename(hp_omim = hp), by = 'identifier') %>%
        na.omit()
      
      list_hpo_genes <- base::split(to_sim_gene$hp_gene, to_sim_gene$gene)      
      list_hpo_diseases <- base::split(to_sim_disease$hp_omim, to_sim_disease$identifier)      
      
      hpo_patient  <- list('patient' = input$chosen_hp)
    
      output_genes <- get_sim_grid(ontology=hpo_dbs, 
                                   term_sets= hpo_patient,
                                   term_sets2 = list_hpo_genes,
                                   term_sim_method = 'resnik') %>%
        t() %>% 
        as_tibble(rownames = 'gene') %>%
        rename(sim_gene = patient) %>%
        mutate(sim_gene = round(sim_gene, 2))

      output_diseases <- get_sim_grid(ontology=hpo_dbs, 
                               term_sets= hpo_patient,
                               term_sets2 = list_hpo_diseases,
                               term_sim_method = 'resnik') %>%
        t() %>% 
        as_tibble(rownames = 'identifier') %>% 
        mutate(identifier = as.double(identifier)) %>%
        rename(sim_mim = patient) %>%
        mutate(sim_mim = round(sim_mim, 2))
    
      return_tbl <- tmp_tbl %>% 
        select(-hp_gene) %>%
        distinct() %>%
        left_join(output_genes, by = 'gene') %>%
        left_join(output_diseases, by = 'identifier') %>%
        na.omit()
        
    }
    
    return_tbl %>%
      left_join(hpo_omim %>% select(identifier, desc, disease_source) %>% distinct(), by = 'identifier') %>%
      select(source, gene, identifier, disease_source, desc, sim_gene, sim_mim ) %>%
      distinct()
    
  })
  
  
  
  
  running_sim_decipher <- reactive({

    
    validate(
      need(length(input$chosen_hp) != 0, "Please, select at least one HP term to describe patient's symptoms")
      
    )
    
    
    tmp_main <- df_overlap_cnvs_running() %>% filter(source == 'decipher') %>% 
      filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>%
      select(id, chrom, start, end, pathogenicity, genotype, variant_class, phenotypes, length_cnv, p_overlap) %>%
      arrange(desc(p_overlap))
    
    validate(
      need(nrow(tmp_main) != 0, "No DECIPHER CNVs found.")
    )
    
    
    tmp_df <- tmp_main %>% 
      filter(phenotypes != '-') %>% 
      select(id, phenotypes) %>%
      separate_rows(phenotypes, sep = '<br>') %>%
      left_join(vector_total_terms %>% select(term, value), by = c('phenotypes' = 'value')) %>%
      na.omit() %>%
      select(id, term)
    
    validate(
      need(nrow(tmp_df) != 0, "No DECIPHER CNVs with phenotype information found.")
    )

    
      list_hpo_decipher <- base::split(tmp_df$term, tmp_df$id)      
      
      hpo_patient  <- list('patient' = input$chosen_hp)
      
      output_decipher <- get_sim_grid(ontology=hpo_dbs, 
                                   term_sets= hpo_patient,
                                   term_sets2 = list_hpo_decipher,
                                   term_sim_method = 'resnik') %>%
        t() %>% 
        as_tibble(rownames = 'id') %>%
        rename(sim_decipher = patient) %>%
        mutate(sim_decipher = round(sim_decipher, 2))

      
      return_tbl <- tmp_main %>% 
        left_join(output_decipher, by = 'id')

      return_tbl
  })
    
    output$dt_running_sim_score <- renderDT({
      
      if (length(input$chosen_hp) != 0) {
      
      length_hp_chosen <- ifelse(length(input$chosen_hp) >= 10, 10, length(input$chosen_hp))

      tmp_df <- running_sim_score()
      tmp_value_total <- p_value_total %>% filter(category != 'cnv' & length_query == length_hp_chosen) %>% select(-length_query)

      
      tmp_df <- tmp_df %>%
        # DISEASE
        left_join(tmp_value_total %>% filter(category == 'disease') %>% mutate(identifier = as.integer(identifier)), by = 'identifier') %>%
        rowwise() %>%
        mutate(tmp_cmp = ifelse(is.na(sim_mim), NA, length(str_split(vector_score, pattern = ', ')[[1]][sim_mim < str_split(vector_score, pattern = ', ')[[1]]]))) %>%
        mutate(is_significant_disease = case_when(
          tmp_cmp < 500 & tmp_cmp > 0 ~ as.character(tmp_cmp/1e4),
          tmp_cmp == 500 ~ '-',
          tmp_cmp == 0 ~ 'lower than 0.0001',
          is.na(tmp_cmp) ~ '-',
        )) %>%
        select(-tmp_cmp, -sim_mim, -vector_score, -category) %>%
        left_join(tmp_value_total %>% filter(category == 'gene') %>% mutate(identifier = as.character(identifier)), by = c('gene' = 'identifier')) %>%
      mutate(tmp_cmp = ifelse(is.na(sim_gene), NA, length(str_split(vector_score, pattern = ', ')[[1]][sim_gene < str_split(vector_score, pattern = ', ')[[1]]]))) %>%
        mutate(is_significant_gene = case_when(
          tmp_cmp < 500 & tmp_cmp > 0 ~ as.character(tmp_cmp/1e4),
          tmp_cmp == 500 ~ '-',
          tmp_cmp == 0 ~ 'lower than 0.0001',
          is.na(tmp_cmp) ~ '-',
        )) %>%
        # select(-tmp_cmp, -sim_gene, -vector_score, -category) %>%
        select(source, gene, disease_source, identifier, desc, is_significant_gene, is_significant_disease)
      
      # tmp_df <- tmp_df %>% replace_na(list(is_significant_gene = '-', is_significant_disease = '-')) %>%
      #   mutate(is_significant_gene = as.numeric(is_significant_gene), 
      #          is_significant_disease = as.numeric(is_significant_disease))

      datatable(tmp_df %>%
                  mutate(identifier = paste0(identifier, '(', disease_source, ')')) %>%
                  select(-disease_source), 
                rownames = FALSE,
                # options = list(dom = 't', pageLength = 25),
                selection = 'single',
                colnames = c('Source', 
                             'Gene', 
                             'Identifier', 
                             'Disease description', 
                             'Similarity score (P-Value) GENE',
                             'Similarity score (P-Value) DISEASE'
                             ))

      } else {
        

      datatable(running_sim_score() %>% replace_na(list(sim_gene = '-', sim_mim = '-')) %>%
                  mutate(identifier = paste0(identifier, '(', disease_source, ')')) %>%
                  select(-disease_source), 
                rownames = FALSE,
                selection = 'single',
                colnames = c('Source', 
                             'Gene', 
                             'Identifier', 
                             'Disease description', 
                             'Similarity score (P-Value) - GENE',
                             'Similarity score (P-Value) - DISEASE'))
        
      }
    })
    
    
    output$decipher_similarity <-  renderDT({
      
      tmp_df <- running_sim_decipher()
      
      length_hp_chosen <- ifelse(length(input$chosen_hp) >= 10, 10, length(input$chosen_hp))
      
      p_value_tmp_df <- p_value_total %>% filter(category == 'cnv' & length_query == length_hp_chosen) %>% 
        select(-length_query, -category)
      
      # test90 <<- tmp_df
      # test91 <<- p_value_tmp_df
      
      tmp_df <- tmp_df %>%
        left_join(p_value_tmp_df, by = c('id' = 'identifier')) %>%
        rowwise() %>%
        mutate(tmp_cmp = ifelse(is.na(sim_decipher), NA, length(str_split(vector_score, pattern = ', ')[[1]][sim_decipher < str_split(vector_score, pattern = ', ')[[1]]]))) %>%
        mutate(is_significant = case_when(
          tmp_cmp < 500 & tmp_cmp > 0 ~ as.character(tmp_cmp/1e4),
          tmp_cmp == 500 ~ '-',
          tmp_cmp == 0 ~ 'lower than 0.0001',
          is.na(tmp_cmp) ~ '-',
        )) %>%
        select(-tmp_cmp, -sim_decipher, -vector_score) %>%
        arrange(desc(is_significant))
      
      
      
      tmp_df <- tmp_df %>% 
        mutate(id = paste0("<a href='", paste0('https://www.ncbi.nlm.nih.gov/clinvar/variation/', id),"' target='_blank'>", id,"</a>")) %>%
        select(p_overlap, id, chrom, start, end, pathogenicity, genotype, variant_class, is_significant, phenotypes)

      
      datatable(tmp_df, 
                extensions = 'Scroller',
                # options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),                escape = FALSE,
                colnames = c('Overlap (%)', 'ID', 'Chrom', 'Start', 'End', 'Pathogenicity', 'Genotype', 'Class', 'Similarity score (P-Value)',
                             'Phenotype'
                             ),
                selection = 'single',
                filter = list(position = 'top'),
                escape = FALSE,
                options = list(
                  autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE,
                  columnDefs = list(list(className = 'dt-center',  targets = c(0:6,8,9, 9)))),  rownames= FALSE)
      
      
      
    })
    

    
  output$n_diseases <- renderUI({
    
    
    number_diseases <- running_sim_score() %>%  pull(identifier) %>%
      unique() %>% length()

    
    tablerInfoCard(
      width = 12,
      value =  paste(number_diseases, 'diseases'),
      status = "primary",
      description = 'Diseases with HPO',
      icon = "database"
    )
    
  })
  
  output$n_genes <- renderUI({
    
    
    n_genes <- running_sim_score() %>%  
      pull(gene) %>%
      unique() %>% 
      length()
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_genes, 'Genes'),
      status = "primary",
      icon = "database",
      description =  'Genes with HPO terms'
      
    )
    
    
    
  })


  output$hpo_assoc_diseases <- renderDT({
    
    validate(
      need(input$dt_running_sim_score_rows_selected != '', "Please, select a row.")
    )
    
    disease_selected <- running_sim_score() %>%
      slice(input$dt_running_sim_score_rows_selected) %>%
      pull(identifier)
    
    
    
    tmp_df <- hpo_omim %>% 
      filter(identifier %in% disease_selected) %>%
      select(identifier, hp, disease_source) %>%
      left_join(hpo_genes %>% select(desc, hp) %>% distinct(), by = 'hp')
    

    tmp_df <- tmp_df %>%  
      mutate(hp = paste0("<a href='", paste0('https://hpo.jax.org/app/browse/term/', hp),"' target='_blank'>", hp,"</a>")) %>%
      
      mutate(identifier = if_else(disease_source == 'OMIM', paste0("<a href='", paste0('https://www.omim.org/entry/', identifier),"' target='_blank'>", identifier,"</a>"),
                                  paste0("<a href='", 
                                         paste0('https://www.orpha.net/consor/cgi-bin/OC_Exp.php?lng=en&Expert=', identifier),"' target='_blank'>", identifier,"</a>"))
      ) %>%
      mutate(identifier = paste0(identifier, ' (', disease_source, ')')) %>%
      select(-disease_source)
                                    
    
    
    datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Disease identifier', 'HPO term', 'Description'))
    
  })
  
  output$hpo_assoc_genes <- renderDT({
    
    validate(
      need(input$dt_running_sim_score_rows_selected != '', "Please, select a row.")
    )
    
    tmp_df <- running_sim_score() %>%
      slice(input$dt_running_sim_score_rows_selected) %>%
      pull(gene)

    tmp_df2 <- hpo_filter() %>% 
      filter(gene %in% tmp_df) %>%
      mutate(hp = paste0("<a href='", paste0('https://hpo.jax.org/app/browse/term/', hp),"' target='_blank'>", hp,"</a>")) %>%
      select(gene, hp, desc) %>%
      distinct()
    
    datatable(tmp_df2, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'HPO term', 'Description'))
    
  })
  
  output$plot_anatomy <- renderHighchart({
    
    
    validate(
        need(length(input$chosen_hp) != 0, "Please, select at least one HP term to describe patient's symptoms"),
        need(input$dt_running_sim_score_rows_selected != '', "Please, select a row.")

    )
    
    
    selected_gene <- running_sim_score() %>% 
      slice(input$dt_running_sim_score_rows_selected) %>% 
      pull(gene)
    
    selected_disease <- running_sim_score() %>% 
      slice(input$dt_running_sim_score_rows_selected) %>% 
      pull(identifier)

    hpo_from_gene <- hpo_genes %>% filter(gene == selected_gene) %>% pull(hp)
    hpo_from_disease <- hpo_omim %>% filter(identifier == selected_disease) %>% pull(hp)
    hpo_from_patient <- input$chosen_hp

    from_gene <- unlist(map(hpo_from_gene, function(x) get_ancestors(hpo_dbs, x))) %>% 
      enframe() %>%
      filter(value %in% anato_df$name) %>%
      count(value) %>% 
      arrange(desc(n)) %>%
      rename(gene = n, name = value) %>%
      mutate(name = as.character(name))
    
    
    from_patient <- unlist(map(hpo_from_patient, function(x) get_ancestors(hpo_dbs, x))) %>% 
      enframe() %>%
      filter(value %in% anato_df$name) %>%
      count(value) %>% 
      arrange(desc(n)) %>%
      rename(Patient = n, name = value) %>%
      mutate(name = as.character(name))
    
    from_disease <-  unlist(map(hpo_from_disease, function(x) get_ancestors(hpo_dbs, x))) %>%
      enframe() %>%
      filter(value %in% anato_df$name) %>%
      count(value) %>%
      arrange(desc(n)) %>%
      rename(disease = n, name = value) %>%
      mutate(name = as.character(name))
    
    plot_df <-  anato_df %>% 
      left_join(from_gene, by = 'name') %>%
      left_join(from_patient,  by = 'name') %>%
      left_join(from_disease, by = 'name') %>%
      replace_na(list("gene" = 0, "Patient" = 0, 'disease' = 0)) %>%
      select(-name) %>%
      pivot_longer(names_to = 'class', values_to = 'valuae', cols = -value) %>%
      mutate(class = if_else(class == 'gene', paste('Gene', paste0('(',selected_gene, ')')), class)) %>%
      mutate(class = if_else(class == 'disease', paste('Disease', paste0('(',selected_disease, ')')), class))

    
    a <- plot_df %>%
      pivot_wider(names_from = class, values_from = valuae ) %>%
        mutate(value = as.factor(value)) %>%
        mutate(value = fct_relevel(value, 'Ear', after = 0 )) %>%
        mutate(value = fct_relevel(value, 'Breast', after = 1 ))
    
    name_columns_tmp <- colnames(a)
    colnames(a) <- c('value','gene', 'patient', 'disease')
    
    highchart() %>% 
      hc_chart(type = "column") %>%
      hc_plotOptions(column = list(stacking = "normal")) %>%
      hc_xAxis(categories = factor(a$value)) %>%
      hc_add_series(name= name_columns_tmp[2],
                    data = a$gene,
                    stack = "Gene") %>%
      hc_add_series(name="Patient",
                    data = a$patient,
                    stack = "Patient") %>%
      hc_add_series(name= name_columns_tmp[4],
                    data = a$disease,
                    stack = "Disease")

  })
  

  # output$plot_similarity_genes <- renderPlot({
  #   
  #   validate(
  #     need(length(input$chosen_hp) > 0, "Please, select at least one HPO term.")
  #   )
  #   
  #   tmp_df <- running_sim_score()
  #   
  # 
  #   if (input$select_sim_gene_disease == 'genes') {
  # 
  #     tmp_df  %>%
  #       select(gene, sim_gene) %>%
  #       arrange(desc(sim_gene)) %>%
  #       distinct() %>%
  #       ggplot(aes(reorder(gene, -sim_gene), sim_gene)) +
  #       geom_col(aes(fill = sim_gene), color = 'black', show.legend = FALSE) +
  #       ylab('Phenotypic similarity score') +
  #       xlab('Genes') +
  #       scale_fill_viridis_c() +
  #       theme_fancy() +
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
  #             axis.title.x = element_text(size = 16),
  #             axis.title.y = element_text(size = 16))
  # 
  #   } else if (input$select_sim_gene_disease == 'diseases') {
  # 
  #     tmp_df  %>%
  #       select(identifier, sim_mim) %>%
  #       distinct() %>%
  #       arrange(desc(sim_mim)) %>%
  #       # filter(similarity_score >= filter_higher_than) %>%
  #       ggplot(aes(reorder(identifier, -sim_mim), sim_mim)) +
  #       geom_col(aes(fill = sim_mim), color = 'black', show.legend = FALSE) +
  #       # coord_flip() +
  #       ylab('Phenotypic similarity score') +
  #       xlab('Disease identifiers') +
  #       scale_fill_viridis_c() +
  #       # scale_x_reverse() +
  #       theme_fancy() +
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #             axis.title.x = element_text(size = 16),
  #             axis.title.y = element_text(size = 16))
  # 
  #   } else if (input$select_sim_gene_disease == 'decipher') {
  #     
  #     
  #     running_sim_decipher()  %>%
  #       select(id, sim_decipher) %>%
  #       distinct() %>%
  #       arrange(desc(sim_decipher)) %>%
  #       # filter(similarity_score >= filter_higher_than) %>%
  #       ggplot(aes(reorder(id, -sim_decipher), sim_decipher)) +
  #       geom_col(aes(fill = sim_decipher), color = 'black', show.legend = FALSE) +
  #       # coord_flip() +
  #       ylab('Phenotypic similarity score') +
  #       xlab('DECIPHER CNVs identifiers') +
  #       scale_fill_viridis_c() +
  #       # scale_x_reverse() +
  #       theme_fancy() +
  #       theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #             axis.title.x = element_text(size = 16),
  #             axis.title.y = element_text(size = 16))
  #     
  #     
  #     
  #     
  #     
  #     
  #   } else {
  #     p1 <- running_sim_decipher() %>% 
  #       filter(variant_class == 'Deletion') %>%
  #       ggplot(aes(p_overlap, sim_decipher, label = id)) + 
  #       geom_point(aes(fill = variant_class), color = 'black', shape = 21, show.legend = FALSE) + 
  #       geom_text_repel() +
  #       theme_minimal() +
  #       labs(title = 'Deletion', x = 'Overlap (%)', y = 'Phenotypic similarity') +
  #       theme(plot.title = element_text(size=25),
  #             axis.title=element_text(size=15,face="bold"))
  #       
  #       p2 <- running_sim_decipher() %>% 
  #         filter(variant_class == 'Duplication') %>%
  #       ggplot(aes(p_overlap, sim_decipher, label = id)) + 
  #         geom_point(aes(fill = variant_class), color = 'black', shape = 21, show.legend = FALSE) + 
  #         geom_text_repel() +
  #         theme_minimal() +
  #         labs(title = 'Duplication', x = 'Overlap (%)', y = 'Phenotypic similarity') +
  #         theme(plot.title = element_text(size=25),
  #               axis.title=element_text(size=15,face="bold"))
  #       
  #       p1 + p2
  #       
  #       
  # 
  #   }
  # 
  # 
  # })
  
  
  output$n_hi<- renderUI({
    
    
    n_hi_yes <- data_selected() %>% filter(hi <= 10) %>% nrow()
    n_total <- nrow(data_selected())  
    
    tablerStatCard(
      value =  paste(n_hi_yes, n_total, sep = '/'),
      title =  paste("Likely to exhibit haploinsufficiency"),
      width = 12
    )

  })

  
  
  output$choose_geno_karyo2 <- renderUI({
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      textInput(
        inputId = "int_end",
        label = "Genomic interval - End",
        value = '153724866')
    }
  })
  
  model_genes_phenotype <- reactive({
 
    chosen_genes <- data_selected() %>% select(gene) %>% pull()

    
   if (length(chosen_genes) == 0) {
     
     mgi_tmp <- tibble()
     
   } else {
    
    mgi_tmp <- mgi %>% filter(gene %in% chosen_genes) %>%
      separate_rows(pheno, sep = ', ') %>%
      rename(term = pheno)
    

    vector_mpo <- mgi_tmp %>% 
      select(term) %>% 
      distinct() %>%
      left_join(mpo_dbs, by = 'term')
    

    mgi_tmp <- mgi_tmp %>% 
      left_join(vector_mpo, by = 'term') %>%
      left_join(data_selected() %>% select(gene, source), by = 'gene') %>%
      select(source, everything())
    
   }

    mgi_tmp
    
  })

  output$agg_model <- renderPlot({
    
    validate(
      need(nrow(model_genes_phenotype()) != 0, '0 genes found.')
    )

    tmp_df <- model_genes_phenotype() %>% 
      count(source, description) %>% arrange(desc(n))

    tmp_df %>%
      ggplot(aes(reorder(description, n), n)) + 
      geom_col(aes(fill = n), color = 'black', show.legend = FALSE) +
      coord_flip() +
      scale_fill_viridis_c() +
      ylab('Number of genes') +
      xlab('Mouse phenotype terms') + 
      theme_fancy() +
      theme(axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 13),
            axis.text.y = element_text(size = 16)) +
      labs(fill = NULL) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      facet_grid(~source)
      

  })

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
      # df_genes <- data_selected()[input$dgenes_rows_all,]
      df_genes <- data_selected()
      
    }
    
    filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
    univ <- hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character()
    
    validate(
      need(length(filtered_genes) != 0, "0 enriched diseases found.")
    )
    
    go_analysis <- enrichGO(gene  = filtered_genes,
                            universe      = univ,
                            OrgDb         = org.Hs.eg.db,
                            ont           = go_chosen,
                            # pAdjustMethod = "BH",
                            pvalueCutoff  = as.numeric(input$sign_vline),
                            readable = TRUE)
    
    # qvalueCutoff  = 0.05)
    
    
    go_analysis <- go_analysis %>% as_tibble() %>% filter(Count > 1)

    
    validate(
      need(nrow(go_analysis) != 0, "0 enriched terms found.")
    )
    
    go_analysis
  })
  
  
  output$func_analysis <- renderPlot({
    
    running_enrich_go() %>%
      mutate(p.adjust = -log10(p.adjust)) %>%
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
      separate_rows(geneID,sep = '/') %>%
      na.omit() %>%
      distinct() %>%
      mutate(pvalue = round(pvalue, 3)) %>%
      mutate(p.adjust = round(p.adjust, 3)) %>%
      mutate(qvalue = round(qvalue, 3)) %>%
      mutate(ID = paste0("<a href='", paste0('http://amigo.geneontology.org/amigo/term/', ID),"' target='_blank'>", ID,"</a>"))
    
    
    datatable(df, escape = FALSE)
    
    
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
      need(ggo %>% as_tibble() %>% select(Count) %>% sum() != 0, "0 terms found. Please reduce the level selected.")
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
    
    
    running_go() %>%
      arrange(desc(Count)) %>%
      filter(Count != 0) %>%
      slice(1:10) %>%
      as_tibble() %>% 
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
    
    
    df <- running_go() %>%
      as_tibble() %>%
      filter(Count != 0) %>%
      arrange(desc(Count)) %>%
      # separate(geneID, sep = '/', into = as.character(1:1000)) %>%
      # gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio) %>%
      # select(-delete) %>%
      separate_rows(geneID, sep ='/') %>%
      na.omit() %>%
      distinct()
    
    datatable(df)
    
  })
  
  running_path <- reactive({
    
    req(isTRUE(input$enable_path_analysis))
    
    # if (is.null(input$dgenes_rows_all)) {
    #   df_genes <- data_selected()
    # } else {
    #   df_genes <- data_selected()[input$dgenes_rows_all,]
    # }
    
    filtered_genes <- data_selected() %>% select(entrez_id) %>% pull()  %>% as.character()
    
    validate(
      need(length(filtered_genes) != 0, "0 genes found.")
    )
    
    if (input$kegg_reactome == 'Reactome') {
      
      pathway_analysis <- enrichPathway(gene=  filtered_genes,
                                        pvalueCutoff = as.numeric(input$sign_vline_path), 
                                        # universe = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
                                        readable= TRUE)  %>% 
        as_tibble()
      # mutate(ID = paste0("<a href='", paste0('https://reactome.org/content/detail/', ID),"' target='_blank'>", ID,"</a>"))
      
    } else {
      pathway_analysis <- clusterProfiler::enrichKEGG(gene= filtered_genes ,
                                                      pvalueCutoff= as.numeric(input$sign_vline_path),
                                                      universe = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character()) %>% 
        as_tibble()
      # mutate(ID = paste0("<a href='", paste0('https://reactome.org/content/detail/', ID),"' target='_blank'>", ID,"</a>"))
    }
    
    pathway_analysis <- pathway_analysis %>%
      filter(Count > 1)
    
    validate(
      need(nrow(pathway_analysis) != 0, "0 enriched pathways found.")
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
  
  
  
  
  
  
  output$df_path  <- renderDT({
    
    df <-  running_path()  %>%
      arrange(desc(Count)) %>%
      # separate(geneID, sep = '/', into = as.character(1:1000)) %>%
      # gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
      # select(-delete) %>%
      separate_rows(geneID, sep = '/') %>%
      distinct()
    
    datatable(df, escape = FALSE)
    
  })

  
  output$func_pathways  <- renderPlot({
    
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
    
    
    # if (is.null(input$dgenes_rows_all)) {
    #   df_genes <- data_selected()
    # } else {
    #   df_genes <- data_selected()[input$dgenes_rows_all,]
    # }
    
    df_genes <- data_selected()

    filtered_genes <- df_genes %>% select(entrez_id) %>% pull()  %>% as.character()
    
    validate(
      need(length(filtered_genes) != 0, "0 enriched diseases found.")
    )
    
    
    
    enrich_dgn <- enrichDO(gene  = filtered_genes,
                           universe      = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
                           pvalueCutoff  = as.numeric(input$pvalue_do),
                           readable = TRUE)
    


    validate(
      need(nrow(enrich_dgn %>% as_tibble() %>% filter(Count > 1)) != 0, "0 enriched terms found.")
    )
    
    enrich_dgn
    
  })
  
  
  output$df_do <- renderDT({

    df <- running_do() %>%
      as_tibble() %>%
      filter(Count > 1) %>%
      arrange(desc(Count)) %>%
      separate_rows(geneID, sep = '/') %>%
      na.omit() %>%
      distinct() %>%
      mutate(pvalue = round(pvalue, 3)) %>%
      mutate(p.adjust = round(p.adjust, 3)) %>%
      mutate(qvalue = round(qvalue, 3)) %>%
      mutate(ID = str_replace(ID, ':', '_')) %>%
      mutate(ID = paste0("<a href='", paste0('https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F', ID),"' target='_blank'>", ID,"</a>"))
    
    
    datatable(df, escape = FALSE)
    
  })
  
  output$func_do  <- renderPlot({
    

    cnetplot(running_do(), foldChange= hgcn_genes$entrez_id, readable = TRUE)
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
  
  
  reading_cnv_file <- reactive({

    req(input$input_geno_karyo == 'Multiple coordinates (NGS)')
    
      validate(
        need(!is.null(input$file_cnv), "Please upload a file.")
      )

    file1 <- input$file_cnv

    data1 <- read_tsv(file1$datapath, col_names = TRUE,
                      col_types = list(chrom = col_character(),
                                       start = col_character(),
                                       end = col_character(),
                                       type = col_character()))
    
    colnames(data1) <- c('chrom', 'start', 'end', 'type')
    
    data1 <- data1 %>% 
      # filter(!across(everything(), is.na)) %>% 
      filter_all(any_vars(!is.na(.))) %>% # remove empty lines (common in txt files)
      mutate(chrom = str_remove(chrom, 'chr')) %>% #remove chr from chrom column
      mutate(start = str_remove_all(start, ',')) %>%
      mutate(end = str_remove_all(end, ',')) %>%
      mutate(start = as.integer(start)) %>%
      mutate(end = as.integer(end)) %>%
      mutate(start = start + 1)  %>% # 0-based to 1-based 
      filter((end - start  + 1) <= 15e6) # remove variants larger than 15 millions b.p
      
      data1
  })
  
  

  
  output$cnv_file <- renderDT({
    
    tmp_tbl <- reading_cnv_file() %>% mutate(length_cnv = end - start + 1)
    
    datatable(tmp_tbl, colnames = c('Chrom', 'Start', 'End','Type', 'Lenght'))
    
  })
  
  cnv_file_to_analyze <- reactive({
    
    tmp_tbl <- reading_cnv_file() %>% select(-type)
    

    if (input$select_all_cnvs == 'yes') {
      
      tmp_tbl
      
    } else {
      
      validate(
        need(length(input$cnv_file_rows_selected) > 0, 'Please, click on the rows to select, at least,
             one variant')
      )
      
      
      tmp_tbl[input$cnv_file_rows_selected,]
    }
    
    
    
    
    
    
  })
  
  df_overlap_cnvs_running <- reactive({

    tmp_df <-  get_perc_overlap(check_cnv_df(), coord_user(), is_patho = TRUE)

    tmp_df
  })

  
  output$df_intersection <- renderDT({
    
    
    df_tmp <- intersection_running() %>%
      select(id, start, end)
    
    datatable(df_tmp, rownames = FALSE, selection = 'single',
              options = list(dom = 't'))
    
  })
  
  df_overlap_cnvs_running_download_patho <- reactive({
    
    
    tmp_df <-  df_overlap_cnvs_running() %>% filter(source == 'decipher')
    
    tmp_df
  })
  
  df_overlap_cnvs_running_download_nonpatho <- reactive({
    
    
    tmp_df <-  df_overlap_cnvs_running() %>% filter(source != 'decipher')
    
    tmp_df
  })

  
  
  output$df_overlap_cnvs <- renderDT({
    
    req(input$select_decipher_clinvar)

    tmp_df <- df_overlap_cnvs_running() %>% filter(source == 'decipher') %>% 
      filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>%
      select(id, chrom, start, end, pathogenicity, genotype, variant_class, phenotypes, length_cnv, p_overlap) %>%
      arrange(desc(p_overlap))
    
    
    validate(
      need(nrow(tmp_df) != 0, "No pathogenic CNVs found.")
    )
    
  
    
    tmp_df <- tmp_df %>%
      rowwise() %>%
      mutate(genes = paste(bed_intersect(hgcn_genes, 
                                         tibble('chrom' = chrom, 'start' = start,
                                                      'end' = end)) %>% pull(gene.x), collapse = ', ')) %>%
      ungroup() %>%
      select(p_overlap, everything())
      # mutate(genes = str_replace_all(genes, ', ', '<br>' ))
    
    
    if (input$select_decipher_clinvar == 'DECIPHER') {

    datatable(tmp_df, 
              extensions = 'Scroller',
              escape = FALSE,
              colnames = c('Overlap (%)', 'ID', 'Chrom', 'Start', 'End', 'Pathogenicity', 'Genotype', 'Class', 'Phenotype',
                           'CNV size', 'Genes overlapping'),
              selection = 'single',
              filter = list(position = 'top'), 
              rownames = FALSE,
              options = list(
                # autoWidth = TRUE,
                deferRender = TRUE, scrollY = 200, scrollX = TRUE, fixedColumns = TRUE))
                # columnDefs = list(list(className = 'dt-center',  targets = c(0:6,8)),
                #                   list(width = '20px', targets = 9))),  rownames= FALSE)
      
    } else {
      
      tmp_df <- running_clinvar_yes_cnv()
      
      validate(
        need(nrow(tmp_df) != 0, "No pathogenic CNVs (ClinVar) found.")
      )
      
      tmp_df <- get_perc_overlap(running_clinvar_yes_cnv(), coord_user(), is_patho = TRUE) %>% 
        select(p_overlap, id, chrom, start, end, variant_class, clinical, disease_identifier, disease_name, gene) %>%
        mutate(disease_identifier = str_replace_all(disease_identifier, 'Human Phenotype Ontology', 'HPO'))



      tmp_df <- tmp_df %>% 
        # mutate(disease_identifier = str_replace_all(disease_identifier, '\\|', '<br>' )) %>%
        # mutate(gene = str_replace_all(gene, ';', '<br>' )) %>%
        # mutate(gene = str_replace_all(gene, ':', '<br>' )) %>%
        mutate(disease_name = if_else(disease_name == 'See cases', '-', disease_name)) %>%
        mutate(id = paste0("<a href='", paste0('https://www.ncbi.nlm.nih.gov/clinvar/variation/', id),"' target='_blank'>", id,"</a>"))
      
      datatable(tmp_df,
                extensions = 'Scroller',
                filter = list(position = 'top'), 
                options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                escape = FALSE, 
                colnames = c('Overlap (%)','ID', 'Chrom','Start','End','Variant Class', 'Clinical significance', 
                             'Disease Identifier', 
                             'Disease name',
                             'Dosage-sensitive genes'), 
                rownames = FALSE
      )
    }
  })
  
  
  output$select_db_no_patho_cnvs <- renderUI({
    
    
    
    vector_tmp <- c(paste('DECIPHER Control',paste0('(', nrow(running_decipher_c()), ')')),
                    paste('DGV',paste0('(', nrow(running_dgv()), ')')),
                    paste('gnomAD v.2.1',paste0('(', nrow(running_gnomad()), ')')))
    
    
    vector_n_dbs <-    split(c('decipher_control', 'dgv', 'gnomad_v2.1'), vector_tmp)
    
    prettyRadioButtons(
      inputId = "select_no_patho_cnv",
      label = '', 
      choices = vector_n_dbs,
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
    
    
  })
  
  
  # output$ui_select_del_dup <- renderUI({
  #   
  #   vector_n_dbs <-    split(c('deletions', 'duplications'), c('Deletions', 'Duplications'))
  #   
  #   prettyRadioButtons(
  #     inputId = "select_del_dup",
  #     label = tags$b("Select one option:"),   
  #     choices =  split(c('deletions', 'duplications'), c('Deletions', 'Duplications')),
  #     inline = TRUE, 
  #     status = "primary",
  #     fill = TRUE
  #   )
  # })
  
  
  
  
  running_dgv <- reactive({
    
    # req(input$start_analysis > 0)

    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'dgv') %>% pull(id)
    
    tmp_df <- dgv_df_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      select(p_overlap, id, chrom, start, end, variantsubtype, reference, pubmedid, method, samplesize, observedgains,
             observedlosses, genes)
    


    
    
    tmp_df
  })
  
  running_decipher_c <- reactive({
    

    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'decipher_control') %>% pull(id)
    
    tmp_df <- decipher_control_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      arrange(desc(p_overlap)) %>%
      select(-source) %>%
      select(p_overlap, id, chrom, start, end, everything())
    
    
    
    
    tmp_df
  })
  
  running_gnomad <- reactive({
    
    # req(input$start_analysis > 0)

    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'gnomad_v2.1') %>% pull(id)
    
    tmp_df <- gnomad_sv_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      select(p_overlap, id, chrom, start, end, svtype, AF)
    
    
    
    tmp_df
  })
  
  output$df_overlap_cnvs_nonpatho <- renderDT({
    
    req(input$select_no_patho_cnv)
    
    
    if (input$select_no_patho_cnv == 'dgv') {
      
      validate(
        need(nrow(running_dgv()) > 0, '0 non-pathogenic CNVs from DGV found.' )
      )

      datatable(running_dgv(), 
                extensions = 'Scroller',
                options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                rownames = FALSE,
                colnames = c('Overlap (%)','ID', 'Chrom', 'Start', 'End', 'Type', 'Reference', 'PMID',
                             'Method', 'Sample size', 'Observed gains', 'Observed losses', 'Genes'
                             ))

    } else if (input$select_no_patho_cnv == 'decipher_control') {
      
      validate(
        need(nrow(running_decipher_c()) > 0, '0 non-pathogenic CNVs from Decipher Population found.' )
      )
      
      
      datatable(running_decipher_c(), 
                extensions = 'Scroller',
                options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                rownames = FALSE,
                colnames = c('Overlap (%)','ID', 'Chrom', 'Start', 'End', 'Deletion Obs.', 'Deletion Freq.',
                             'Deletion (SE)', 'Duplication Obs.', 'Duplication Freq.', 'Duplication (SE)',
                             'Observations', 'Frequency', 'Standard Error', 'Type','Sample size', 'Study'
                             ))
    } else {
      
      validate(
        need(nrow(running_gnomad()) > 0, '0 non-pathogenic CNVs from gnomAD v.2.1 found.')
      )
      
      datatable(running_gnomad(), 
                extensions = 'Scroller',
                options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                rownames = FALSE,
                colnames = c('Overlap (%)', 'ID', 'Chrom', 'Start', 'End', 'Type', 'Allele Frequency'))
    }
    
    
  })
  
  upload_raw <- reactive({
    
    req(input$upload_bed_file)
    
    x <- input$upload_bed_file

    tmp_path <-  x$datapath
    tmp_df <- read_tsv(tmp_path, col_names = c('chrom', 'start', 'end'))
    
    
  })
  
  
  running_upload <- reactive({
    
    req(upload_raw())
    
    tmp_df <- upload_raw() %>% mutate(size = end - start) %>%
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
  
  
  # 
  # output$df_upload <- renderDT({
  #   
  #   req(running_upload())
  # 
  #   tmp_df <- running_upload() %>%
  #     mutate(score = pmap(list(chrom, start, end), function(a,b,c) get_model_score(a,b,c, hgcn_genes, model1))) %>%
  #     rowwise() %>%
  #     mutate(n_genes = score[[2]],
  #            score = score[[1]])
  #   
  # 
  #   
  #   datatable(tmp_df, colnames = c('Chromosome', 'Start', 'End', 'Size', 'Score', 'nº genes'), rownames = FALSE)
  #   
  #   
  # })
  running_gwas <- reactive({
    
    # req(input$start_analysis > 0)
    

    tmp_df <- gwas_variants %>%
      # rename(chrom = CHR_ID, start = CHR_POS) %>%
      # mutate(end = start) %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      rename(pos = start) %>%
      select(-end) %>%
      mutate(pubmed_id = str_extract(LINK, '\\d{8}')) %>%
      select(-LINK) %>%
      select(-contains('source')) %>%
      distinct()

  })
  

  
  
  running_de_novo <- reactive({
    
    tmp_df <- denovo %>%
      rename(start = Position) %>%
      mutate(end = start) %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      rename(position = start) %>%
      select(-end) %>%
      distinct()

      
    tmp_df
    
  })
  
  running_clinvar <- reactive({
    
    tmp_df <- clinvar_variants %>%
      mutate(chrom = as.character(chrom)) %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      select(-contains('source')) %>%
      distinct() %>%
      mutate(disease_identifier = str_replace_all(disease_identifier, ',', '<br>')) %>%
      select(id, everything())

    tmp_df
    
  })
  
  running_clinvar_no_cnv <- reactive({
    
    running_clinvar() %>% filter(length_cnv < 50) 
  })
  
  running_clinvar_yes_cnv <- reactive({
    
    running_clinvar() %>% filter(length_cnv >= 50) 
  })
  
  
  output$ui_select_clinvar_gwas <- renderUI({
    
    req(running_clinvar_no_cnv())
    req(running_gwas())
    
    
    vector_tmp <- c(paste('ClinVar',paste0('(', nrow(running_clinvar_no_cnv()), ')')),
                    paste('GWAS', paste0('(', nrow(running_gwas()), ')')))
    
    
    vector_n_dbs <-    split(c('clinvar', 'gwas'), vector_tmp)
    
    prettyRadioButtons(
      inputId = "select_clinvar_gwas",
      label = '', 
      choices =  vector_n_dbs,
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
    
  })
  
  output$df_variants <- renderDT({
    
    req(input$select_clinvar_gwas)
    
    if (input$select_clinvar_gwas == 'clinvar') {
      
      validate(
        need(nrow(running_clinvar_no_cnv()) != 0, '0 ClinVar variants found.')
      )
    
      tmp_df <- running_clinvar_no_cnv() %>% 
        # mutate(disease_identifier = str_replace_all(disease_identifier, '\\|', '<br>' )) %>%
        select(id, variant_class, chrom, start, end, reference, alternative, clinical, gene, disease_name, disease_identifier) %>%
        mutate(id = paste0("<a href='", paste0('https://www.ncbi.nlm.nih.gov/clinvar/variation/', id),"' target='_blank'>", id,"</a>")) %>%
        mutate(disease_name = str_replace_all(disease_name, '\\|', '<br>' )) %>%
        mutate(disease_name = str_replace_all(disease_name, ';', '<br>' )) %>%
        mutate(variant_class = if_else(variant_class == 'single nucleotide variant', 'SNV', variant_class))
      

      datatable(tmp_df, 
                extensions = 'Scroller',
                options = list(autoWidth = TRUE, scrollY = 200, scroller = TRUE, scrollX = TRUE, fixedColumns = TRUE),
                filter = list(position = 'top'), 
                escape = FALSE, 
                colnames = c('Clinvar ID','Class', 'Chrom','Start','End', 'Ref', 'Alt','Clinical significance', 'Gene',
                                                     'Disease name', 
                                                     'Disease identifier'), 
                rownames = FALSE
      )
      
    } else {
      
      validate(
        need(nrow(running_gwas()) != 0, '0 GWAS variants found.')
      )
      
      tmp_df <- running_gwas() %>%
        mutate(pubmed_id = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pubmed_id),"' target='_blank'>", pubmed_id,"</a>")) %>%
        mutate(SNPS = paste0("<a href='", paste0('https://www.ebi.ac.uk/gwas/variants/', SNPS),"' target='_blank'>", SNPS,"</a>"))
        
      datatable(tmp_df, 
                extension = 'Scroller',
                options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE),
                filter = list(position = 'top'), 
                escape = FALSE,
                colnames = c('Variant','Chrom', 'Position','Intergenic', 'Reported gene', 'Disease or trait', 'Link study'),
                rownames = FALSE)

    }

  })
  
  
  output$df_de_novo <- renderDT({
    
    validate(
      need(nrow(running_de_novo()) != 0, '0 de novo variants found.')
    )
    
    tmp_df <- running_de_novo() %>% 
      mutate(CaddScore = replace_na(CaddScore, '-'),
             LofScore = replace_na(LofScore, '-')) %>%
      select(-contains('source')) %>%
      mutate(PubmedID = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', PubmedID),"' target='_blank'>", PubmedID,"</a>"))

    datatable(tmp_df, 
              extensions = 'Scroller',
              options = list(deferRender = TRUE, scrollY = 200, scroller = TRUE),
              escape = FALSE, 
              colnames = c('Chrom', 'Position','Gene', 'Phenotype', 'Study name', 'PubmedID', 'Function Class', 
                                                   'CADD score', 'LoF score'), rownames = FALSE)
    
    
  })

  output$ui_intersect <- renderUI({
    
    
    req(input$df_intersection_rows_selected)
    
    actionBttn(
      inputId = "take_intersect",
      label = "Load this region!",
      color = "warning",
      style = "material-flat",
      size = 'sm',
      block = TRUE
    )
    
    
    
    
  })
  
  
  temporal_network <- reactive({
    
    validate(
      need(nrow(data_selected()) != 0, "No protein-coding genes found.")
    )
    
    
    df_nodes <- data_selected() %>% 
      select(gene, source) %>% 
      mutate(id =  row_number() - 1)
    
    df_links <- interactions_db %>% 
      filter(from %in% df_nodes$gene | to %in% df_nodes$gene) %>%
      left_join(df_nodes %>% select(-source), by = c( 'from' = 'gene')) %>% 
      rename(id_from = id) %>%
      left_join(df_nodes %>% select(-source), by = c( 'to' = 'gene')) %>% 
      rename(id_to = id) %>% 
      na.omit() %>%
      select(id_from, id_to)
    
    result <- list('nodes' = df_nodes, 'links' = df_links)
    result
    
  })
  
  output$output_select_gene <- renderUI({
    
    genes_selected <- setNames(as.list(temporal_network()[[1]]$gene), temporal_network()[[1]]$gene)
    
    
               selectInput(
                 "network_ppi_select_gene", '', genes_selected
               )
    
  })
  
  
  output$network_ppi <- renderForceNetwork({

    
    validate(
      need(nrow(temporal_network()[[2]]) != 0, "No protein-protein interactions found.")
    )
    
    ColourScale <- 'd3.scaleOrdinal()
            .domain(["CNV", "Enhancer", "lncRNAs", "miRNAs", "TFs", "TADs"])
           .range(["#66C2A5","#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#ae7373"]);'
    
    
    if (input$filter_by_gene_ppi == 'No') {
      
      forceNetwork(Links = as.data.frame(temporal_network()[[2]]), Nodes = as.data.frame(temporal_network()[[1]]),
                   Source = "id_from", Target = "id_to",
                   NodeID = "gene", 
                   zoom = TRUE,
                   Group = 'source',
                   opacity = 1,
                   opacityNoHover = 1,
                   colourScale = JS(ColourScale))
      
    } else {

      
      req(input$network_ppi_select_gene)
      
      tmp_node <-  temporal_network()[[1]] %>% filter(gene == input$network_ppi_select_gene) %>% pull(id)
      
      gene_chosen2 <-  temporal_network()[[2]] %>% filter(id_from == tmp_node)
      
      
      validate(
        need(nrow(gene_chosen2) != 0, paste("No protein-protein interactions with the ", input$network_ppi_select_gene, 'gene'))
      )
      

      gene_chosen1 <-  temporal_network()[[1]] %>% 
        filter(id %in% c(gene_chosen2 %>% pull(id_from), gene_chosen2 %>% pull(id_to))) %>%
        mutate(id2 = row_number() - 1)

      
      
      gene_chosen2 <- gene_chosen2 %>% 
        left_join(gene_chosen1, by = c('id_from' = 'id')) %>% 
        rename(id_from_real = id2) %>% 
        left_join(gene_chosen1, by = c('id_to' = 'id')) %>% 
        rename(id_to_real = id2)
      
      gene_chosen1 <- gene_chosen1 %>% select(-id) %>% rename(id = id2)

      forceNetwork(Links = as.data.frame(gene_chosen2), Nodes = as.data.frame(gene_chosen1),
                   Source = "id_from_real", Target = "id_to_real",
                   NodeID = "gene", 
                   zoom = TRUE,
                   Group = 'source',
                   opacity = 1,
                   colourScale = JS(ColourScale))
 
    }
    

    
    
  })
  

  
  
  output$frequency_network <- renderPlot({
    
    validate(
      need(nrow(temporal_network()[[2]]) != 0, "No protein-protein interactions found.")
    )

    selected_nodes <- temporal_network()[[2]] %>% count(id_from)
    
    temporal_network()[[1]] %>% 
      left_join(selected_nodes, by = c('id' = 'id_from')) %>% 
      na.omit() %>% 
      arrange(desc(n)) %>%
      slice(1:10) %>%
      mutate(gene = paste0(gene, ' (', source, ')')) %>%
      ggplot(aes(reorder(gene, -n), n)) +
        geom_col(aes(fill =n), color = 'black', show.legend = FALSE) +
        scale_color_viridis_c() +
        scale_y_continuous(breaks = pretty_breaks()) +
        theme_minimal() +
      labs(y = 'Number of interactions', x = 'Gene')

  })
  
  output$drugbank_df <- renderDT({
    
    vector_filt_genes <- data_selected() %>% pull(gene)
    
    
    
    
    tmp_df <- drugbank %>%
      filter(gene %in% vector_filt_genes) %>%
      mutate(uniprot_id = paste0("<a href='", 
                                 paste0('https://www.uniprot.org/uniprot/', uniprot_id),"' target='_blank'>", uniprot_id,"</a>")) %>%
      mutate(drug_id = paste0("<a href='", 
                              paste0('https://www.drugbank.ca/drugs/', drug_id),"' target='_blank'>", drug_id,"</a>"))
    
    
    validate(
      need(nrow(tmp_df) > 0, "0 approved drug targets found.")
    )
    
    datatable(tmp_df, escape = FALSE,
              colnames = c('Target gene', 'Description', 'Uniprot ID', 'DrugBank ID', 'Name', 'Synonyms'))
  })
  
  
  output$doc_chosen <- renderUI({
    
    
    if (input$select_doc_element == 'overview') {
      
      includeMarkdown('../doc/documentation.Rmd')
      
      
    } else if (input$select_doc_element == 'browser') {
      
      includeMarkdown('../doc/browser_compatibility.Rmd')
      
    } else if (input$select_doc_element == 'contact') {
    
      includeMarkdown('../doc/contact.Rmd')
    
    } else if (input$select_doc_element == 'tutorials') {
      
      includeMarkdown('../doc/tutorials.Rmd')
      
    } else if (input$select_doc_element == 'faqs') {
      
      includeMarkdown('../doc/faqs.Rmd')
      
    } else if (input$select_doc_element == 'versions') {
      
      includeMarkdown('../doc/versions.Rmd')
      
    } else if (input$select_doc_element == 'installation') {
      
      includeMarkdown('../doc/installation.Rmd')
    } else if (input$select_doc_element == 'terms_of_use') {
    
    includeMarkdown('../doc/terms_of_use.Rmd')
  }
  })
  
  
  output$ui_tad <- renderUI({
    
    
    req(input$df_tads_rows_selected)
    
    actionBttn(
      inputId = "include_tad_genes",
      label = "Load this region!",
      color = "warning",
      style = "material-flat",
      size = 'sm',
      block = TRUE
    )
    
    
    
    
  })
  
    # 
    # observeEvent(input$link_to_docu, {
    #   tags$a(href = "#shiny-tab-docu")
    # })
  
  # observeEvent(input$link_to_docu, {
  #   newvalue <- "docu"
  #   updateTabItems(session, "panels", newvalue)
  # })
  
  down_file1 <- reactive({
    
    
    tibble('chrom' = c('1','16', '6', '20', 'X', '22'),
           'start' = c(1002999, 34813719, 34813719, 1499999, 1309999, 1199999),
           'end' = c(1009000, 36278623, 36278630, 1600000, 1370000, 1800000),
           'type' = c('DEL', 'DEL', 'DEL', 'DUP', 'DUP', 'DEL'))
    
  })


  output$download_file_1 <- downloadHandler(
    filename = 'file_example.tsv',
    content = function(file) {
      
      write_tsv(down_file1(), file)
    }
  )
  
  
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