function(input, output, session) {
  
  
  
  
  # res_auth <- secure_server(
  #   check_credentials = check_credentials(credentials)
  # )
  # 
  # output$auth_output <- renderPrint({
  #   reactiveValuesToList(res_auth)
  # })
  
  #azteca
  # waitress <- waitress$new(theme = "overlay-percent") # call the waitress
  # 
  # observeEvent(input$start_analysis, {
  #     waitress$
  #       start()$
  #       auto(percent = 5, ms = 150) # increase by 5 percent every 150 milliseconds
  #     Sys.sleep(3.5)
  #     waitress$hide()
  #   })
  # 
  
  # check_popup <- reactive({
  #   
  #   start_coordinates <- coord_user()[1]
  #   end_coordinates <- coord_user()[2]
  #   
  #   test02399 <<- start_coordinates
  #   test0231 <<- end_coordinates
  #   
  #   
  #   test_length <- if((end_coordinates - start_coordinates + 1) < 0)
  #   
  #   
  # observeEvent(((coord_user()[2] - coord_user()[1] + 1) < 0), {
  #   shinyalert("The end of the genomic interval is lower than the start", type = "error")
  # })
  
  # })
  
  
  
  
  
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
    shinyjs::reset('mirnas_selected_tfs')
    
    #
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
  
  observeEvent(input$reset_pheno_analysis, {
    
    
    shinyjs::reset('chosen_hp')
  })
  
  map_blacklist <- reactive({

    blacklist_encode %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap)
    
    
    
  })
  
  
  output$check_blacklist <- renderUI({
    
    tmp_df <- map_blacklist() 
    
    
    req(nrow(tmp_df) > 0)
    
    if (tmp_df %>% count(class) %>% nrow() > 1) {
      
      say_this <- 'High signal and low mappability region(s) overlapping with the CNV'
      
    } else if (str_detect(tmp_df$class[1], 'High')) {
      
      say_this <- 'High signal region(s) overlapping with the CNV'
      
    } else {
      say_this <- 'Low mappability region(s) overlapping with the CNV'
    }
    
    
    tablerInfoCard(
      width = 12,
      value =  nrow(tmp_df),
      status = "danger",
      icon = "database",
      description =  say_this)
    
    
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
    
    
    
    coord_start <- as.numeric(input$int_start)
    coord_end <-  as.numeric(input$int_end)
    coord_chrom <- input$input_chrom
    
    
    # observeEvent(TRUE, {
    #   # Show a modal when the button is pressed
    #   shinyalert("Oops!", "Something went wrong.", type = "error")
    # })
    
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      
      tbl_output <- tibble('chrom' = coord_chrom, 'start' = coord_start,
                           'end' = coord_end)
      
      
    } else if (input$input_geno_karyo == 'G banding') {
      
      df_tmp <- chromPlot::hg_cytoBandIdeo %>%
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

      test24124120419241491924912941299412 <<- cnv_file_to_analyze()

      tbl_output <- cnv_file_to_analyze()
      
      
    }
    
    
    test2020 <<- tbl_output
    
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
      select(-startdelete, -enddelete, -.overlap)

    
    df_output
    
  })
  
  # n_hpo_unique <- reactive({
  #   
  #   req(input$start_analysis > 0)
  #   
  # 
  #   # if (is.null(input$dgenes_rows_all)) {
  #   #   df_genes <- data_selected()
  #   # } else {
  #   #   df_genes <- data_selected()[input$dgenes_rows_all,]
  #   # }
  #   
  #   
  #  symbol_genes <-  data_selected_prev() %>% select(gene) %>% pull()
  #   hpo_genes_filter <- hpo_genes %>% 
  #     filter(gene %in% symbol_genes)
  #   
  #   test9999999 <<- hpo_genes_filter
  # 
  # })
  # 
  check_hp_genes <- reactive({
    
    req(input$start_analysis > 0)
    
    # go <- rols::Ontology("hp")
  
    # if (is.null(input$dgenes_rows_all)) {
    #   df_genes <- data_selected()
    # } else {
    #   df_genes <- data_selected()[input$dgenes_rows_all,]
    # }
    
    
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
  
  output$check_genes_hp <- renderUI({
    
    req(input$start_analysis > 0)
    
    df_genes <- data_selected()
    
    
    
    tablerStatCard(
      value =  nrow(check_hp_genes()),
      title = "Number of genes associated with the phenotype(s)",
      # trend = -10,
      width = 12
    )
    
  })
  
  
  # test_hpo <- reactive({
  #   
  #   
  #   req(input$start_analysis > 0)
  #   # req(input$run_pheno_analysis)
  #   
  #   
  #   genes_selected <- data_selected() %>% select(gene) %>% pull()
  #   
  #   hpo_patient  <- list('patient' = input$chosen_hp)
  #   
  #   tmp_df <- hpo_genes %>% 
  #     filter(gene %in% genes_selected) %>%
  #     select(gene, hp)
  #   
  #   
  #   if (length(input$chosen_hp) == 0) {
  #     
  #     
  #     mat_test <- tmp_df %>%
  #       select(gene) %>%
  #       distinct() %>%
  #       mutate(similarity_score = NA)
  #     
  #     
  #     
  #   } else {
  #     
  # 
  #     
  #     to_p_value <- tmp_df %>% count(gene, name = 'n_freq')
  #     
  #     genes_sample <- base::split(tmp_df$hp, tmp_df$gene)
  #     
  #     
  #     
  #     mat_test <- get_sim_grid(ontology=hpo_dbs, 
  #                              term_sets= hpo_patient,
  #                              term_sets2 = genes_sample,
  #                              term_sim_method = 'resnik') %>%
  #       t()
  #     
  #     
  #     colnames(mat_test) <- c('value')
  #     
  #     mat_test %>% as_tibble(rownames = 'gene') %>% 
  #       mutate(value = round(value, 2)) %>%
  #       rename(similarity_score = value)
  #     
  #   }
  #   
  # 
  #   
  # })
  
  

  calculate_sim_disease <- reactive({
    
    req(input$start_analysis > 0)
    
    
    
    tmp <- data_selected() %>% select(gene) %>% pull()
    
    input_mim_disease <- omim %>% 
      filter(gene %in% tmp) %>%
      pull(MIM_pheno_number) %>%
      unique()
    
    disease_genes_selected <- hpo_omim %>% 
      filter(mim_disease %in% input_mim_disease) %>%
      select(mim_disease, hp) 
    
    if (input$input_inheritance != 'Any') {
      
      selected_mim <- disease_genes_selected %>% 
        filter(hp %in% input$input_inheritance ) %>% 
        pull(mim_disease)
      
      validate(
        need(length(selected_mim) != 0, "No diseases found with the mode of inheritance selected.")
      )
      
      test0000000 <<- selected_mim
      disease_genes_selected <- disease_genes_selected %>%
        filter(mim_disease %in% selected_mim)
      test007 <<- disease_genes_selected
    }
    
    disease_genes_selected <- base::split(disease_genes_selected$hp, disease_genes_selected$mim_disease)
    
    
    
    if (length(input$chosen_hp) == 0) {
      

      mat_test <- hpo_omim %>% 
        filter(mim_disease %in% input_mim_disease) %>%
        select(mim_disease) %>%
        distinct() %>%
        mutate(similarity_score = NA)
      
      

    } else {
    
    hpo_selected  <- list(input$chosen_hp)
    
    
    
    mat_test <- get_sim_grid(ontology=hpo_dbs, 
                             term_sets= hpo_selected,
                             term_sets2 = disease_genes_selected,
                             term_sim_method = 'resnik') %>%
      t()  %>% 
      as_tibble(rownames = 'mim_disease') %>%
      mutate(mim_disease = as.numeric(mim_disease)) %>%
      mutate(V1 = round(V1, 3)) %>%
      rename(similarity_score = V1)
    
    }
    
    mat_test
    
    
  })
  
  output$df_check_hp_genes <- renderDataTable({
    
    req(input$start_analysis > 0)
    

    validate(
      need(nrow(check_hp_genes()) != 0, "Not genes found.")
    )
    
    datatable(check_hp_genes(), options = list(
      pageLength = 5))
  })
  
  # output$plot_upset <- renderPlot({
  #   
  #   test771 <<- check_hp_genes()
  # get_upset(check_hp_genes())
  #   
  #   
  # })
  
  output$n_hp_chosen <- renderUI({
    
    req(input$start_analysis > 0)
    
    
    tablerInfoCard(
      width = 12,
      value =  paste(length(input$chosen_hp), 'HPO terms'),
      status = "primary",
      icon = "database",
      description =  'Select more terms in the panel below'
      
    )
    
    # tablerStatCard(
    #   value =  length(input$chosen_hp),
    #   title = "Clinical features selected",
    #   # trend = -10,
    #   width = 12
    # )
    
  })
  
  output$n_cnv_patho <- renderUI({
    
    req(input$start_analysis > 0)
    
    data_tmp <- check_cnv_df() %>% 
      filter(source == 'decipher' & pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>% 
      nrow()
    
    tablerStatCard(
      value =  data_tmp,
      title = "Pathogenic CNVs (DECIPHER)",
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  data_tmp,
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of articles found in Pubmed associated with deletions",
    #   width = 12
    # )
    
    
  })
  
  output$n_cnv_nopatho <- renderUI({
    
    req(input$start_analysis > 0)
    
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
    
    
    if ((input$input_geno_karyo == 'Multiple coordinates' | input$input_geno_karyo == 'Genomic coordinates')) {

      tmp_query <-  chromPlot::hg_cytoBandIdeo %>%
        rename(chrom = Chrom, start = Start, end = End) %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        mutate(chrom_name = paste0(chrom, Name)) %>%
        select(chrom, Name, chrom_name) %>%
        mutate(result = paste('(','(', 'chromosome', chrom,'AND', '(', Name ,'OR', chrom_name,  ')', ')')) %>%
        pull(result) %>%
        paste(collapse = ' OR ') %>%
        paste('AND deletion AND homo sapiens')
        

      
      query_region <- tmp_query
      

    } else {
      
      chrom_tmp <- paste('chromosome',  coord_user() %>%  pull(chrom))
      band_tmp <-  input$input_karyotype
      band2_tmp <- paste0( coord_user() %>%  pull(chrom), band_tmp)
      
      query_region <- paste(chrom_tmp,'AND','(', band_tmp,'OR', band2_tmp,')', 'AND deletion AND homo sapiens')
      
    }
    
    query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 200 )

  })
  
  
  
  # query_pubmed_dup <- reactive({
  #   
  #   req(input$start_analysis > 0)
  #   
  #   start_coordinates <- coord_user()[1]
  #   end_coordinates <- coord_user()[2]
  #   chrom_coordinates <- coord_user()[3]
  #   
  #   
  #   if (input$input_geno_karyo == 'Genomic coordinates') {
  # 
  #     query_region <-  chromPlot::hg_cytoBandIdeo %>%
  #       filter(Chrom %in% chrom_coordinates) %>%
  #       mutate(keep = map2_chr(Start, End, function(x,y) c(start_coordinates, end_coordinates) %overlaps% c(x,y))) %>%
  #       filter(keep == TRUE) %>%
  #       select(Name) %>%
  #       pull() %>%
  #       map_chr(function(x) paste0(chrom_coordinates, x)) %>%
  #       paste0(collapse = ' OR ') %>%
  #       paste('AND duplication', sep = ' ')
  #     
  #   } else {
  #     query_region <- paste0(chrom_coordinates, input$input_karyotype)
  #     
  #   }
  #   
  #   test67 <<- query_region
  #   
  #   query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 200 )
  #   
  #   
  #   
  # })
  
  output$n_pubmed_del <- renderUI({
    
    # tablerInfoCard(
    #   value =  length(query_pubmed_del()[['ids']]),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of articles found in Pubmed associated with deletions",
    #   width = 12
    # )
    
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
  
  
  
  # output$omim_assoc <- renderDataTable({
  # 
  #   ids_query <- c(query_pubmed_dup()[['ids']], query_pubmed_del()[['ids']])
  #   test0123456 <<- ids_query
  # 
  #   ids_query <- test0123456
  #   query_link <- entrez_link(db= 'omim', id= ids_query, dbfrom="pubmed", by_id = TRUE)
  #   query_link <- query_link$links[['pubmed_omim_calculated']] %>% as.numeric()
  # 
  #   validate(
  #     need(query_link != 0, 'No OMIM entries')
  #   )
  # 
  # 
  #   query_tmp <- entrez_summary(db="omim", id= query_link)
  # 
  #   title <- unname(map_chr(query_tmp, function(x) x[["title"]]))
  # 
  # 
  # 
  # 
  #   tmp_df <- tibble('title' = title,'id_omim' = query_link) %>%
  #     mutate(id_omim = paste0("<a href='", paste0('https://www.omim.org/entry/', id_omim),"' target='_blank'>", id_omim,"</a>"))
  # 
  #   datatable(tmp_df, escape = FALSE, colnames = c('Title', 'ID OMIM'))
  # 
  # 
  # })
  
  
  
  
  omim_assoc <- reactive({
    
    
    if (input$select_del_dup == 'deletions') {
      
      ids_query <- query_pubmed_del()[['ids']]
    } else {
      ids_query <- query_pubmed_dup()[['ids']]
    }
    # 
    # 
    # ids_query <- c(query_pubmed_dup()[['ids']], query_pubmed_del()[['ids']])
    # 
    # ids_query <- ids_query[1:200]
    
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
    
    test12311 <<- tmp_df 

    validate(
      need(nrow(tmp_df  %>% filter(omim_assoc != '')) != 0, 'No articles associated with OMIM entries.')
    )
    
    
    tmp_df
  })
  
  
  
  output$gene_assoc <- renderDataTable({
    
    ids_query <- c(query_pubmed_del()[['ids']], query_pubmed_dup()[['ids']])
    
    test01234 <<- ids_query
    query_link <- entrez_link(db= 'gene', id= ids_query, dbfrom="pubmed")
    test4141 <<- query_link
    query_link <- query_link$links[['pubmed_gene']] %>% as.numeric()
    
    validate(
      need(query_link != 0, 'No gene entries')
    )
    
    test1311111 <<- query_link
    query_tmp <- entrez_summary(db="gene", id= test1311111)
    
    title <- unname(map_chr(query_tmp, function(x) x[["name"]]))
    id_gene <- unname(map_chr(query_tmp, function(x) x[["uid"]]))
    id_chrom <-unname(map_chr(query_tmp, function(x) x[["chromosome"]]))
    
    organism <- unname(map_chr(query_tmp, function(x) x[['organism']][['scientificname']]))
    
    
    tmp_df <- tibble('title' = title, 'organism' = organism, 'chromosome' = id_chrom, 'id_gene' = id_gene) %>%
      mutate(id_gene = paste0("<a href='", paste0('https://www.ncbi.nlm.nih.gov/gene/', id_gene),"' target='_blank'>", id_gene,"</a>")) 
    
    datatable(tmp_df, escape = FALSE, colnames = c('Title', 'Organism', 'Chromosome', 'ID GENE'), filter = 'top')
    
    
  })
  
  
  # output$medgen_assoc <- renderDataTable({
  #   
  #   ids_query <- c(query_pubmed_del()[['ids']], query_pubmed_dup()[['ids']])
  #   query_link <- entrez_link(db= 'medgen', id= ids_query, dbfrom="pubmed")
  #   test4141 <<- query_link
  #   query_link <- query_link$links[['pubmed_medgen']] %>% as.numeric()
  #   
  #   validate(
  #     need(FALSE, 'TBA')
  #   )
  #   
  #   test1311111 <<- query_link
  #   query_tmp <- entrez_summary(db="medgen", id= test1311111)
  #   
  #   title <- unname(map_chr(query_tmp, function(x) x[["name"]]))
  #   id_medgen <- unname(map_chr(query_tmp, function(x) x[["uid"]]))
  # 
  #   
  #   tmp_df <- tibble('title' = title, 'id_gene' = id_medgen)
  #   
  #   datatable(tmp_df, escape = FALSE, colnames = c('Title', 'id'))
  #   
  #   
  # })
  
  query_pubmed_dup <- reactive({
    
    
    
    
    if ((input$input_geno_karyo == 'Multiple coordinates' | input$input_geno_karyo == 'Genomic coordinates')) {
      
      tmp_query <-  chromPlot::hg_cytoBandIdeo %>%
        rename(chrom = Chrom, start = Start, end = End) %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        mutate(chrom_name = paste0(chrom, Name)) %>%
        select(chrom, Name, chrom_name) %>%
        mutate(result = paste('(','(', 'chromosome', chrom,'AND', '(', Name ,'OR', chrom_name,  ')', ')')) %>%
        pull(result) %>%
        paste(collapse = ' OR ') %>%
        paste('AND duplication AND homo sapiens')
      
      
      test91241412 <<- tmp_query
      
      # if (length(tmp_query) > 1) tmp_query %>% coll 
      query_region <- tmp_query
      
      
    } else {
      
      chrom_tmp <- paste('chromosome',  coord_user() %>%  pull(chrom))
      band_tmp <-  input$input_karyotype
      band2_tmp <- paste0( coord_user() %>%  pull(chrom), band_tmp)
      
      query_region <- paste(chrom_tmp,'AND','(', band_tmp,'OR', band2_tmp,')', 'AND duplication AND homo sapiens')
      
    }
    
    # test91214 <<- query_region
    
    query_pubmed <- entrez_search(db="pubmed", term= query_region, retmax = 200 )

  })
  
  output$n_pubmed_dup <- renderUI({
    
    
    tablerStatCard(
      value =   length(query_pubmed_dup()[['ids']]),
      # status = "primary",
      # icon = 'book',
      title = "Articles found in Pubmed associated with duplications",
      width = 12
    )
    
    
    # tablerInfoCard(
    #   value =   length(query_pubmed_dup()[['ids']]),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of articles found in Pubmed associated with duplications",
    #   width = 12
    # )
    
  })
  
  
  output$n_mortality <- renderUI({
    
    tmp_n <- model_genes_phenotype() %>%  filter(description == 'mortality/aging') %>% nrow()
    

    tablerStatCard(
      value =  tmp_n,
      title = "Genes associated with mortality/aging phenotype",
      width = 12
    )
    
    
  })
  
  output$n_embryo <- renderUI({
    
    tmp_n <- model_genes_phenotype() %>%  filter(description == 'embryo') %>% nrow()

    tablerStatCard(
      value =   tmp_n,
      title = "Genes associated with embryonic phenotype",
      width = 12
    )
    
  })
  
  
  # run_arules <- reactive({
  #   
  #   # test913 <<- as.numeric(input$type_cnv)
  #   # test914 <<- as.numeric(input$denovo_yes_no)
  #   
  #   tmp_variant <- as.numeric(input$type_cnv) # Deletion 1 - Duplication 0
  #   tmp_inheritance <- as.numeric(input$denovo_yes_no) # De novo consitutive 1 - Other 0
  #   
  #   # tmp_start <- as.numeric(coord_user()[1])
  #   # tmp_end <- as.numeric(coord_user()[2])
  #   # tmp_chrom <- coord_user()[3]
  #   
  #   # tmp_start <- 34813719
  #   # tmp_end <- 36278623
  #   # tmp_chrom <- '4'
  #   
  #   tmp_out  <- check_app_cnv( 1, 1, tmp_variant, tmp_inheritance,
  #                              tmp_chrom, tmp_start, tmp_end)  %>% 
  #     select(-id, -clinical) %>%
  #     rename(maximum_pli =  max_pli)
  #   
  #   test99912 <<- tmp_out
  #   object <- model_cba
  #   newdata <- tmp_out
  #   type <- 'class'
  #   method <- model_cba$method
  #   
  #   
  #   
  #   newdata <- arules::discretizeDF(tmp_out, lapply(object$discretization, 
  #                                                   FUN = function(x) list(method = "fixed", breaks = x)))
  #   
  #   newdata <- as(newdata, "transactions")
  #   newdata <- arules::recode(newdata, match = lhs(object$rules))
  #   rulesMatchLHS <- is.subset(lhs(object$rules), newdata, sparse = (length(newdata) * 
  #                                                                      length(rules(object)) > 150000))
  #   dimnames(rulesMatchLHS) <- list(NULL, NULL)
  #   
  #   matched_rules <- rulesMatchLHS %>% as_tibble(name = 'V1') %>% 
  #     mutate(id = row_number()) %>% 
  #     filter(V1 == TRUE) %>% 
  #     pull(id)
  #   
  #   validate(
  #     need( length(matched_rules) != 0, "No association rules matching")
  #   )
  #   
  #   test900410410 <<- matched_rules
  #   
  #   tmp_output <- bind_cols(quality(model_cba$rules[matched_rules]), 
  #                           labels(model_cba$rules[matched_rules],  itemSep = " + ") %>% 
  #                             enframe(name = 'rule')) %>%
  #     slice(1)
  #   
  #   conf_result <- tmp_output$confidence
  #   
  #   tmp_object <- tmp_output$value %>% str_split('=>')
  #   
  #   test0202 <<- tmp_object
  #   
  #   prediction_result <- tmp_object[[1]][2] %>% 
  #     str_remove('\\{clinical=') %>%
  #     str_remove('\\}')
  #   
  #   rules_result <- tmp_object[[1]][1] %>% 
  #     str_split('\\+')
  #   
  #   test0101 <<- list('prediction' = prediction_result, 
  #                     'rules' = rules_result)
  #   
  #   output_list <- list('prediction' = prediction_result, 
  #                       'rules' = rules_result,
  #                       'confidence' = round(conf_result,1)*100)
  #   
  # })
  
  
  
  # output$rules_arules <- renderDT({
  #   
  # 
  #   df_rules <- run_arules()[[2]][[1]] %>% 
  #     enframe(name = 'rule') %>%
  #     separate(value, sep = '=', into = c('def', 'value')) %>%
  #     mutate(def = str_remove(def, '\\{'),
  #            def = str_replace(def, 'n_cnv_syndromes', 'Nº of CNV syndromes overlapping'),
  #            def = str_replace(def, 'disease_variants', 'Nº of disease variants overlapping'),
  #            def = str_replace(def,  'embryo_mouse', 'Nº of genes associated with lethality in embryonic mouse'),
  #            def = str_replace(def, 'n_genes_hpo', 'Nº of genes associated with HPO terms'),
  #            def = str_replace(def, 'n_genes', 'Nº of genes'),
  #            def = str_replace(def,  'type_variant', 'Type of CNV'),
  #            def = str_replace(def,  'type_inheritance', 'Type of inheritance'),
  #            def = str_replace(def, 'nonpatho_cnv', 'Nº of non-pathogenic CNVs overlapping'),
  #            def = str_replace(def, 'patho_cnv', 'Nº of pathogenic CNVs overlapping'),
  #            def = str_replace(def, 'disease_genes', 'Nº disease genes'),
  #            def = str_replace(def, 'maximum_pli', 'Maximum pLI score found'),
  #            def = str_replace(def, 'n_tf', 'Nº of transcription factors (TFs)')
  #     ) %>%
  #     mutate(value = str_replace(value, '-Inf', '0'),
  #            value = str_replace(value, '\\[', 'from '),
  #            value = str_replace(value, ',', ' to '),
  #            value = str_replace(value, '\\)', ''),
  #            value = str_replace(value, '\\]', ''),
  #            value = str_replace(value, '\\}', ''),
  #            value = str_replace(value, 'from 0.5 to  Inf', 'higher than 0'),
  #            value = str_replace(value, 'from 0 to 0.5', 'equal to 0')
  #            ) %>%
  #     mutate(rule = paste0("<b>", rule, "</b>"),
  #            def = paste0("<b>", def, "</b>"))
  #   
  #   datatable(df_rules, 
  #             escape = FALSE,
  #             rownames = FALSE, 
  #             colnames = c('Rule', 'Description', 'Value'),
  #             options = list(dom = 't',
  #                            columnDefs = list(list(className = 'dt-center', targets = c(0:2)))), 
  #             class = 'cell-border stripe')
  #   
  #   
  # })
  
  output$ui_run_arules  <- renderUI({
    
    
    
    tablerStatCard(
      value =   run_arules()[[1]],
      title = "Prediction - association rules",
      width = 12
    )
    
  })
  
  output$ui_confidence  <- renderUI({
    
    tablerStatCard(
      value =   paste0(run_arules()[[3]], '%'),
      title = paste("of CNVs of DECIPHER that follow these rules are considered:",
                    run_arules()[[1]]),
      width = 12
    )
    
  })
  
  
  running_pubmed_del <- reactive({
    
    
    
    # validate(
    #   need(length(query_pubmed()) != 0, "Please, select a gene in the datatable.")
    # )
    
    # test21 <<- query_pubmed()
    query_tmp <- entrez_summary(db="pubmed", id= query_pubmed_del()[['ids']])
    
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
    
    df_output <- df_output %>% 
      mutate(n_cites = as.integer(ifelse(n_cites == '', 0, n_cites)))
    
    df_output
    
    
  })
  
  running_pubmed_dup <- reactive({
    
    req(query_pubmed_dup())
    
    # validate(
    #   need(length(query_pubmed()) != 0, "Please, select a gene in the datatable.")
    # )
    
    # test21 <<- query_pubmed()
    query_tmp <- entrez_summary(db="pubmed", id= query_pubmed_dup()[['ids']])
    
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
  
  output$del_dup_pubmed <- renderDataTable({
    
    req(running_pubmed_del())
    req(input$select_del_dup)
    
    
    if (input$select_del_dup == 'deletions') {
      
      tmp_df <- running_pubmed_del()
      
    } else {
      
      tmp_df <- running_pubmed_dup()
      
    }
    
    test010 <<- tmp_df
    
    
    if(input$only_omim == 'Yes') {
      
      
      
      tmp_df <- tmp_df %>% 
        select(pmid, everything()) %>%
        left_join(omim_assoc(), by = c('pmid' = 'pubmed_id')) %>%
        filter(omim_assoc != '') %>%
        mutate(omim_assoc = paste0("<a href='", paste0('https://www.omim.org/entry/', omim_assoc),"' target='_blank'>", omim_assoc,"</a>")) %>%
        # mutate(tmp_col = 'tmp_cols') %>%
        # pivot_wider(id_cols = pmid, values_from = omim_assoc, names_from = tmp_col) %>%
        mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>"))
      
      
      
      
      
      
      vector_colnames <- c('PMID', 'Title','First author', 'Last author', 'N°cites','Journal', 'Published date', 'OMIM entries')
      
    } else {
      
      tmp_df <- tmp_df %>% 
        select(pmid, everything()) %>%
        mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>"))
      
      
      vector_colnames <- c('PMID', 'Title','First author', 'Last author', 'N°cites','Journal', 'Published date')
      
    }
    

    
    datatable(tmp_df, rownames = FALSE, filter = 'top', selection = 'single', escape = FALSE,
              
              colnames = vector_colnames,
              options = list(
                pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                selection = 'single'
                # columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
    
  })
  
  output$abstract_del_pubmed <- renderDataTable({
    
    validate(
      need(input$del_dup_pubmed_rows_selected != '', "Please, select an article from the table.")
    )
    
    # test931 <<- input$dup_pubmed_rows_selected
    # paper_selected <- running_pubmed_del() %>% slice(input$del_pubmed_rows_selected)
    
    if (input$select_del_dup == 'deletions') {
      
      paper_selected <- running_pubmed_del() %>% slice(input$del_dup_pubmed_rows_selected)
      
    } else {
      
      paper_selected <- running_pubmed_dup() %>% slice(input$del_dup_pubmed_rows_selected)
      
    }
    
    
    fetch.pubmed <- entrez_fetch(db = "pubmed", id = paper_selected %>% pull(pmid), rettype = "xml", parsed = T)
    # Extract the Abstracts for the respective IDS.  
    abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
      xmlValue(xmlChildren(x)$Abstract))
    
    tmp_df <- tibble(Title = paper_selected$title, Abstract = abstracts[[1]]) %>% gather('Category', 'Info')
    tmp888 <<- tmp_df
    
    datatable(tmp_df, rownames = FALSE, colnames = '',
              options = list(dom = 't'))
    
  })
  
  
  output$dup_pubmed <- renderDataTable({
    
    
    tmp_df <- running_pubmed_dup()
    tmp_df <- tmp_df %>% mutate(pmid = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pmid),"' target='_blank'>", pmid,"</a>")) 
    
    
    
    
    
    datatable(tmp_df, rownames = FALSE, filter = 'top', selection = 'single', escape = FALSE,
              
              colnames = c('Title','First author', 'Last author', 'N°cites','Journal', 'Published date', 'PMID' ),
              options = list(
                pageLength = 5, autoWidth = TRUE, style = 'bootstrap', list(searchHighlight = TRUE),
                selection = 'single'
                # columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
    
  })
  
  
  
  output$abstract_pubmed <- renderDataTable({
    
    validate(
      need(input$dup_pubmed_rows_selected != '', "Please, select an article from the datatable.")
    )
    
    paper_selected <- running_pubmed_dup() %>% slice(input$dup_pubmed_rows_selected)
    
    
    fetch.pubmed <- entrez_fetch(db = "pubmed", id = test45 %>% pull(pmid), rettype = "xml", parsed = T)
    
    abstracts = xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
      xmlValue(xmlChildren(x)$Abstract))
    
    tmp_df <- tibble(title = paper_selected$title, abstract = abstracts[[1]]) %>% gather('Category', 'Info')
    
    datatable(tmp_df, rownames = FALSE, escape = FALSE,
              options = list(dom = 't'))
    
  })
  
  
  
  data_selected_prev <- reactive({
    
    
    req(coord_user())
    

    data_raw <- hgcn_genes %>% mutate_if(is.factor, as.character)
      
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
        
        # Genes  mapped in CNV
        genes_cnv <- data_selected_prev()
        genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # Genes NOT mapped in CNV
        if (is.null(input$df_enhancer_rows_all)) {
          enhancers_df <- prev_enhancer()
        } else {
          # enhancers_df <- prev_enhancer()[input$df_enhancer_rows_all,]
          enhancers_df <- prev_enhancer()
          #
        }
        genes_no_cnv <- enhancers_df %>% select(gene) %>% distinct() %>% pull()
        genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
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
  
  data_selected_tads <- reactive({
    
    if (!is.null(input$tads_on_off)) {
      if (input$tads_on_off) {
        
        # Genes  mapped in CNV
        genes_cnv <- data_selected_prev()
        genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # Genes NOT mapped in CNV
        if (is.null(input$df_tads_rows_all)) {
          tads_df <- prev_tads()
        } else {
          tads_df <- prev_tads()[input$df_tads_rows_all,]
        }
        genes_no_cnv <- tads_df %>% select(genes_not_cnv) %>% pull() %>% str_split(pattern = ', ') %>% unlist() %>% unique()
        genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        # ADAPT IT WHEN ADDING OMIM OR OTHERS!!!
        table_output <- hgcn_genes %>% filter(gene %in% genes_no_cnv) %>%
          select(-vg, -ensembl_gene_id) %>%
          mutate(source = 'TAD')
        
        
      } else {
        table_output <- tibble()
      }
    } else {
      table_output <- tibble()
    }
    
  })
  
  data_selected_mirnas <- reactive({
    
    
    
    if (!is.null(input$mirnas_on_off)) {
      if (input$mirnas_on_off) {
        
        # Genes  mapped in CNV
        genes_cnv <- data_selected_prev()
        genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # Genes NOT mapped in CNV
        
        if (is.null(input$df_mirna_rows_all)) {
          mirnas_df <- mirna_raw()
        } else {
          mirnas_df <- mirna_raw()
          # mirnas_df <- mirna_raw()[input$df_mirna_rows_all,]
        }
        
        
        genes_no_cnv <- mirnas_df %>% select(gene_symbol) %>% distinct() %>% pull()
        genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        
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
  
  data_selected_tfs <- reactive({
    
    
    
    if (!is.null(input$tfs_on_off)) {
      if (input$tfs_on_off) {
        
        # Genes  mapped in CNV
        genes_cnv <- data_selected_prev()
        genes_cnv <- genes_cnv %>% select(gene) %>% pull()
        # Genes NOT mapped in CNV
        
        if (is.null(input$df_enhancer_rows_all)) {
          tfs_df <- tf_raw()
        } else {
          # tfs_df <- tf_raw()[input$tf_df_rows_all,]
          tfs_df <- tf_raw()
          
        }
        
        
        genes_no_cnv <- tfs_df %>% select(target) %>% distinct() %>% pull()
        genes_no_cnv <- genes_no_cnv[! genes_no_cnv %in% genes_cnv]
        
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
  
  
  
  data_selected <- reactive({
    
    # observeEvent({input$start_analysis,
    #   
    #   data_selected_prev()
    #   
    #   
    # })
    
    
    table_output <- data_selected_prev() %>% 
      mutate(source = 'CNV') %>% 
      bind_rows(data_selected_enhancers()) %>%
      bind_rows(data_selected_tads()) %>%
      bind_rows(data_selected_tfs()) %>%
      bind_rows(data_selected_mirnas()) %>%
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
  
  
  running_cnv_syndromes <- reactive({
    

    tmp_df <- syndromes_total %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      replace_na(replace = list(variant_class = '-', phenotypes = '-')) %>%
      select(chrom, start, end, syndrome_name, variant_class, phenotypes, source) %>%
      get_perc_overlap((coord_user()))
    
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
  
  output$cnv_syndromes <- renderDataTable({
    
    req(running_cnv_syndromes())
    
    validate(
      need(nrow(running_cnv_syndromes()) != 0, "No CNV syndromes found.")
    )
    
    
    if (input$select_cnv_syndrome == 'decipher') {
      
      tmp_df <- running_cnv_syndromes() %>% 
        filter(source == 'decipher') %>% 
        select(chrom, start, end, syndrome_name, variant_class, phenotypes, p_overlap)
      
      
      datatable(tmp_df, rownames = FALSE,
                colnames = c('Chrom', 'Start', 'End', 'CNV syndrome name', 'Variant class', 'Phenotypes', 'Overlap (%)'))
      
    } else {
      
      tmp_df <- running_cnv_syndromes() %>% filter(source == 'clingen') %>% select(chrom, start, end, syndrome_name, p_overlap, -source)

      
      datatable(tmp_df, rownames = FALSE,
                colnames = c('Chrom', 'Start', 'End', 'CNV syndrome name', 'Overlap (%)'))
      
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
    
    test63321311321312312441241266 <<- uspset_df
    uspset_df
  })
  
  
  output$plot_upset_disease <- renderPlot({
    
    test0001 <<- running_upset_disease()
    
    get_upset(test0001, gene = TRUE)
    
    
    
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
    
    test63321311321312312441241266 <<- uspset_df
    uspset_df
  })
  
  output$plot_upset_disease_reg <- renderPlot({
    
    test0001 <<- running_upset_disease_reg()
    
    get_upset(test0001, gene = TRUE)
    
    
    
  })
  
  
  
  output$dgenes <- renderDataTable({
    
    server <- TRUE
    
    
    
    data_input <- data_selected()  %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      select(-source) %>%
      select(band, gene, disease, orphanet, dev, clingen, omim, gwas, p_overlap) %>%
      filter(disease == 'Yes')
    # mutate(n_evidences = sample(1:4, n(), replace = TRUE)) %>%
    # select(band, gene, n_evidences)
    
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
  
  output$dgenes_reg <- renderDataTable({
    
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
    
    
    test <- split(options_sources, 
                  toupper(options_sources) %>% str_replace('_', ' ') %>% str_replace('DEV', 'DECIPHER') %>%
                    str_replace('GENOMICS ENGLAND', 'GENOMICS ENGLAND PanelApp'))
    
    pickerInput(
      inputId = "select_source",
      label = "", 
      choices = test
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
    
    
    test <- split(options_sources, 
                  toupper(options_sources) %>% str_replace('_', ' ') %>% str_replace('DEV', 'DECIPHER') %>%
                    str_replace('GENOMICS ENGLAND', 'GENOMICS ENGLAND PanelApp'))
    
    pickerInput(
      inputId = "select_source_reg",
      label = "", 
      choices = test
    )
  })
  
  
  output$select_gene_disease <- renderDataTable({
    
    req(input$dgenes_rows_selected)
    
    validate(
      need(input$dgenes_rows_selected != '', 'Please, select a disease gene on the left panel.')
    )
    
    
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
        select(gene, Name6, Name, OrphaNumber, SourceOfValidation ) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.orpha.net/consor/cgi-bin/Disease_Genes.php?lng=EN&data_id=16132&Disease_Disease_Genes_diseaseGroup=', gene),"' target='_blank'>", gene,"</a>")) %>%
        mutate(OrphaNumber = paste0("<a href='", 
                                    paste0('https://www.orpha.net/consor/cgi-bin/Disease_Search.php?lng=EN&data_id=8648&Disease_Disease_Search_diseaseGroup=', OrphaNumber),"' target='_blank'>", OrphaNumber,"</a>") )
      
      
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene',
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
  
  output$select_gene_disease_reg <- renderDataTable({
    
    req(input$dgenes_reg_rows_selected)
    
    validate(
      need(input$dgenes_reg_rows_selected != '', 'Please, select a disease gene on the left panel.')
    )
    
    
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
        select(gene, Name6, Name, OrphaNumber, SourceOfValidation ) %>%
        mutate(gene = paste0("<a href='", 
                             paste0('https://www.orpha.net/consor/cgi-bin/Disease_Genes.php?lng=EN&data_id=16132&Disease_Disease_Genes_diseaseGroup=', gene),"' target='_blank'>", gene,"</a>")) %>%
        mutate(OrphaNumber = paste0("<a href='", 
                                    paste0('https://www.orpha.net/consor/cgi-bin/Disease_Search.php?lng=EN&data_id=8648&Disease_Disease_Search_diseaseGroup=', OrphaNumber),"' target='_blank'>", OrphaNumber,"</a>") )
      
      
      
      datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Gene',
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
  
  output$dgenes_no_disease <- renderDataTable({
    
    server <- TRUE
    
    tmp_df <-  data_selected() 
    
    data_input <- tmp_df %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV') %>%
      select(-source) %>%
      select(band, gene, disease, essent, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
             ncgerp, p_overlap) %>%
      filter(disease == 'No')
    
    validate(
      need(nrow(data_input) > 0, "0 no disease genes.")
    )
    
    
    tmp_output <- datatable(data_input, rownames = FALSE, 
                            colnames = c('Band', 'Gene', 'Disease', 'Essential', 'pLI', 'RVIS', 'CCR', 'HI', 'GDI', 'SnIPRE', 'ncRVIS',
                                         'ncGERP', 'Overlap(%)'),
                            filter = list(position = 'top'), 
                            selection = 'single'
                            # options = list(dom = 't')
    ) %>%
      formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp', 'essent'), color = styleInterval(94, c('weight', '#ff7f7f'))) %>%
      formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f')))
    # formatStyle(c('disease', 'haplo', 'triplo', 'omim', 'dev', 'fda', 'gwas'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
    
    
    tmp_output
    
    
    
  })
  
  
  
  
  # outputOptions(output, "dgenes", suspendWhenHidden= TRUE)
  
  output$choose_reg_region <- renderUI({
    
    
    tmp_data <- data_selected() %>% 
      filter(source != 'CNV') %>%
      filter(disease != 'Yes') %>%
      count(source) %>%
      mutate(source = as.character(source)) %>%
      na.omit()
    
    
    test99999987 <<- tmp_data
    
    
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
  
  
  
  output$genes_from_reg_regions <- renderDataTable({
    
    # if (!is.null(input$enhancers_on_off)) {
    #   if (input$enhancers_on_off) {
    
    
    validate(
      need(!is.null(input$select_reg_region), "0 non-disease target genes selected.")
    )
    
    server <- TRUE
    data_input <- data_selected() %>% filter(source ==  input$select_reg_region) %>%
      select(-start, -end, -chrom) %>%
      select(band, gene, disease, essent, pLI, rvis, ccr, hi, gdi, snipre, ncrvis, 
             ncgerp) %>%
      filter(disease == 'No')
    
    datatable(data_input, rownames = FALSE, filter = list(position = 'top'),
              colnames = c('Band', 'Gene', 'Disease', 'Essential', 'pLI', 'RVIS', 'CCR', 'HI', 'GDI', 'SnIPRE', 'ncRVIS',
                           'ncGERP')) %>%
      formatStyle(c('pLI', 'rvis', 'hi', 'gdi', 'snipre', 'ncrvis', 'ncgerp'), color = styleInterval(94, c('weight', '#ff7f7f'))) %>%
      formatStyle(c('ccr'), color = styleInterval(1, c('weight', '#ff7f7f'))) %>%
      formatStyle(c('disease'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
    
    # formatStyle(c('disease', 'haplo', 'triplo', 'omim', 'dev', 'fda', 'gwas'), color = styleEqual(c('No', 'Yes'), c('weight', '#ff7f7f')))
    
    # 
    #   }
    # }
    # )
  })
  
  # output$genes_from_tads <- renderDataTable({
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
    
    req(input$gtex_gene_tissue)
    
    
    if (input$gtex_gene_tissue == 'Gene' ) {
      
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
    
    gs<-GeneSet(geneIds= vector_genes,
                organism="Homo Sapiens",
                geneIdType=SymbolIdentifier())
    
    output <- teEnrichment(inputGenes = gs, 
                           rnaSeqDataset = if_else(input$tissue_expression_dbs == 'gtex', 1, 3))
    
    
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
  
  output$tissue_hpa <- renderDT({
    
    #   req(input$input_gene_tissue)
    
    vector_genes <- data_selected() %>% pull(gene)
    tmp_df <- hpa %>% filter(gene %in% vector_genes)
    
    if (input$gene_yes_no == 'Yes') {
      req(input$input_gene_tissue)
      filtered_gene <- input$input_gene_tissue
      test99 <<- filtered_gene
      tmp_df <- tmp_df %>% filter(gene == !!filtered_gene)
    }
    
    if (input$tissue_yes_no == 'Yes') {
      req(input$input_tissue)
      filtered_tissue <- input$input_tissue
      tmp_df <- tmp_df %>% filter(tissue == !!filtered_tissue)
    }
    
    # 
    # observe({
    #   if (input$choice == 'Hello') {
    #     getStatus <- 'Hi there'
    #   }
    
    datatable(tmp_df, filter = 'top')
    
    
    
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
    
    df_tmp <- model_genes_phenotype() %>% 
      select(-entrez_id, -gene_mouse) %>%
      mutate(mgi = paste0("<a href='", paste0('http://www.informatics.jax.org/marker/', mgi),"' target='_blank'>", mgi,"</a>")) %>%
      mutate(mpo_id = paste0("<a href='", paste0('http://www.informatics.jax.org/vocab/mp_ontology/', mpo_id),"' target='_blank'>", mpo_id,"</a>"))
    
    
    datatable(df_tmp, escape = FALSE, rownames = FALSE, 
              colnames = c('Human ortholog gene', 'Mouse gene', 'Mouse phenotype id', 'Phenotype description'),
              options = list(
                columnDefs = list(list(className = 'dt-center', targets = 0:3) )))
    
    
    
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
  
  # output$dot_comparison <- renderPlotly({
  #   
  #   validate(
  #     need(FALSE, "TBA")
  #   )
  #   
  #   validate(
  #     need(input$dgenes_rows_selected != '', "Please, select a gene in the datatable.")
  #   )
  #   
  #   name_gene_filtered <- data_selected() %>% slice(input$dgenes_rows_selected) %>% select(gene) %>% pull()
  #   
  #   p <- data_selected() %>%
  #     
  #     ggplot(aes(pLI, rvis)) +
  #     geom_point(col = "steelblue") +
  #     theme_fancy() +
  #     gghighlight(gene == name_gene_filtered)
  #   
  #   p <- ggplotly(p)
  #   p <- plotly_build(p)
  #   
  #   p$x$data[1][[1]]$text <- paste0('Gene symbol: ', data_selected() %>% pull(gene), 
  #                                   "<br>",
  #                                   'pLI: ', data_selected()  %>% pull(pLI),
  #                                   "<br>",
  #                                   'RVIS: ', data_selected()  %>% pull(rvis)
  #                                   
  #   )
  #   
  #   p$x$data[[2]]$text <- paste0('Gene symbol: ', data_selected()  %>% slice(input$dgenes_rows_selected) %>% pull(gene), 
  #                                "<br>",
  #                                'pLI: ', data_selected() %>% slice(input$dgenes_rows_selected) %>% pull(pLI),
  #                                "<br>",
  #                                'RVIS: ', data_selected() %>% slice(input$dgenes_rows_selected)  %>% pull(rvis)
  #   )
  #   
  #   test1898 <<- p
  #   p
  #   
  #   
  # })
  
  
  
  
  
  
  output$data <- renderTable({
    mtcars[, c("mpg", input$variable), drop = FALSE]
  }, rownames = TRUE)
  
  
  
  plot_chrom_react <- reactive({
    
  req(input$input_geno_karyo != 'Multiple coordinates')

    start_coordinates <- coord_user() %>% pull(start)
    end_coordinates <- coord_user() %>% pull(end)
    chrom_coordinates <- coord_user() %>% pull(chrom)
    chrom_coordinates <- paste0('chr', chrom_coordinates)
    
    
    ideoTrack <- IdeogramTrack(genome="hg19", chromosome= chrom_coordinates)
    plotTracks(ideoTrack, from= start_coordinates , to = end_coordinates, showBandId=TRUE,
               cex.bands=0.5)
    
    # plotKaryotype(chromosomes = input_chr, plot.type = 2) %>%
    # kpDataBackground(data.panel = 1)
    # kpAddBaseNumbers() %>%
    # kpRect(chr= input_chr, x0= data_input$start_position, x1= data_input$end_position, y0=0.2, y1=0.4)
    
    
  })
  
  
  # fisher_running <- reactive({
  #   
  #   req(input$start_analysis > 0)
  #   
  #   
  #   list_panel_names <- panel_total %>% select(Level4) %>% distinct() %>% pull()
  #   input_genes <- data_selected() %>% select(gene) %>% pull()
  #   
  #   fisher_result <- tibble()
  #   
  #   for (i in 1:length(list_panel_names)) {
  #     # print(i)
  #     input_test <- matrix(rep(NA, 4), ncol = 2, dimnames = list(c('in_panel', 'no_panel'), c('col1', 'col2') ))
  #     
  #     input_test[1,][2] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% pull(gene) %in% input_genes %>% sum()
  #     input_test[1,][1] <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% nrow() - input_test[1,][2] 
  #     input_test[2,][2] <- length(input_genes) - input_test[1,][2] 
  #     input_test[2,][1] <- 19146 - input_test[2,][2] - input_test[1,][1] - input_test[1,][2]
  #     
  #     name_genes <- panel_total %>% filter(Level4 == !!list_panel_names[i]) %>% pull(gene)
  #     name_genes <- name_genes[name_genes %in% input_genes]
  #     if (length(name_genes) == 0) {
  #       name_genes <- '-'
  #     }
  #     
  #     tmp_tibble <- tibble(name_panel = list_panel_names[i], 
  #                          p_value = fisher.test(input_test, alternative = 'greater')$p.value,
  #                          genes  = name_genes,
  #                          gene_ratio = paste(as.character(input_test[1,][2]), '/',
  #                                             as.character(input_test[1,][2] + input_test[1,][1]))
  #     )
  #     fisher_result <- rbind(fisher_result, tmp_tibble)
  #   }
  #   
  #   fisher_result %>% arrange(p_value) %>% mutate(p_value = round(p_value, 4))
  # 
  # })
  
  
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
  # 
  # fisher_running_two <- reactive({
  #   
  #   req(input$start_analysis > 0)
  #   
  #   
  #   name_dbs <- c('omim', 'clinvar', 'haplo', 'triplo', 'dev', 'fda', 'gwas', 'essent')
  #   input_genes <- data_selected() %>% select(gene) %>% pull()
  #   # input_genes <- test2019 %>% select(gene) %>% pull()
  #   fisher_result <- tibble()
  #   
  #   for (i in 1:length(name_dbs)) {
  #     # i <- 8
  #     print(i)
  #     input_test <- matrix(rep(NA, 4), ncol = 2, dimnames = list(c('in_panel', 'no_panel'), c('col1', 'col2') ))
  #     
  #     input_test[1,][2] <- hgcn_genes %>% filter(get(name_dbs[i]) == 'Yes') %>% pull(gene) %in% input_genes %>% sum()
  #     input_test[1,][1] <- hgcn_genes %>% filter(get(name_dbs[i]) == 'Yes') %>% nrow() - input_test[1,][2] 
  #     input_test[2,][2] <- length(input_genes) - input_test[1,][2] 
  #     input_test[2,][1] <- 19146 - input_test[2,][2] - input_test[1,][1] - input_test[1,][2]
  #     test312321 <<- input_test
  #     tmp_tibble <- tibble(name_panel = name_dbs[i], 
  #                          p_value = fisher.test(input_test, alternative = 'greater')$p.value,
  #                          gene_ratio = paste(as.character(input_test[1,][2]), '/',
  #                                             as.character(input_test[1,][2] + input_test[1,][1]))
  #     )
  #     test44444 <<- tmp_tibble
  #     fisher_result <- rbind(fisher_result, tmp_tibble)
  #   }
  #   
  #   fisher_result %>% arrange(p_value) %>% mutate(p_value = round(p_value, 4))
  #   
  # })
  
  # output$df_fisher_two <- renderDataTable({
  #   
  #   datatable(fisher_running_two(), colnames = c('Database name', 'p.value', 'Gene ratio'))
  #   
  # })
  
  
  
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
    
    # test1000 <<-coord_user()
    
    plot_chrom_react()
    
    
  })
  
  output$choose_geno_karyo1 <- renderUI({
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      numericInput(
        inputId = "int_start",
        label = "Genomic interval - Start",
        value = 34813719)
      
    } else {
      
      karyotype_filtered <- as.list(chromPlot::hg_cytoBandIdeo %>% filter(Chrom == input$input_chrom) %>% select(Name))
      
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
    
    # test1111111 <- data_selected()
    # if (is.null(input$dgenes_rows_all)) {
    #   final_number <- test1111111
    # } else {
    #   final_number <- test1111111[input$dgenes_rows_all,]
    # }
    
    
    
    
    # tablerStatCard(
    #   value = nrow(data_selected()),
    #   title = HTML(paste0("Number of genes<br/>", name_region)),
    #   # trend = 19192,
    #   width = 12
    # )
    
    tablerInfoCard(
      width = 12,
      value =  paste(nrow(data_selected()), 'Genes'),
      status = "primary",
      icon = "database",
      description = 'Number of genes - CNV')
    # description = HTML(paste0("Number of genes<br/>", name_region))
    
    
    
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
  
  
  
  
  
  # output$ref_user_filter_genes <- renderUI({
  # 
  #  req(input$dgenes_rows_all)  
  #  req(length(input$dgenes_rows_all) != nrow(data_selected() %>% 
  #                                              select(-start_position, -end_position, -chrom) %>%
  #                                              filter(source == 'CNV')))
  # 
  #   test18 <<- length(input$dgenes_rows_all) 
  #   test20 <<-  nrow(data_selected() %>% 
  #                      select(-start_position, -end_position, -chrom) %>%
  #                      filter(source == 'CNV'))
  #   
  # 
  #   if (input$input_geno_karyo == 'Genomic coordinates') {
  # 
  #     name_region <- paste0('chr',coord_user()[3], ':', input$int_start, '-', input$int_end)
  # 
  #   } else {
  #     name_region <- paste0(coord_user()[3], input$input_karyotype)
  # 
  #   }
  # 
  #   tablerInfoCard(
  #     width = 12,
  #     value = paste0(length(input$dgenes_rows_all), '/', nrow(data_selected()), " genes"),
  #     status = "warning",
  #     icon = "crop",
  #     description =  'Filtered genes'
  #   )
  # 
  # 
  # })
  
  # output$ref_user_filter_genes2 <- renderUI({
  #   
  #   req(length(input$dgenes_rows_all) != nrow(data_selected()))
  #   
  #   
  #   if (input$input_geno_karyo == 'Genomic coordinates') {
  #     
  #     name_region <- paste0('chr',coord_user()[3], ':', input$int_start, '-', input$int_end)
  #     
  #   } else {
  #     name_region <- paste0(coord_user()[3], input$input_karyotype)
  #     
  #   }
  #   
  #   tablerInfoCard(
  #     width = 12,
  #     value = paste0(length(input$dgenes_rows_all), '/', nrow(data_selected()), " genes"),
  #     status = "warning",
  #     icon = "crop",
  #     description =  'Filtered genes'
  #   )
  #   
  #   
  # })
  # 
  
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
      status = "primary",
      icon = "database",
      # description =  name_region
      description = 'Total number of genes'
    )
    
    
  })
  
  # outputOptions(output, "dgenes", suspendWhenHidden = FALSE)
  
  # observeEvent(input$start_analysis, {
  # 
  #   dt_test_proxy <- dataTableProxy("dgenes", session = shiny::getDefaultReactiveDomain(),
  #                                   deferUntilFlush = TRUE)
  #   
  #   replaceData(dt_test_proxy, tmp_output)
  # 
  # })
  
  
  output$ref_user_genes_cnv <- renderUI({
    
    # shinyjs::reset('dgenes_rows_all')
    
    
    rows_selected <- input$dgenes_rows_all
    
    
    # observeEvent(input$start_analysis, {
    #   
    #   # shinyjs::hide("dgenes")
    #   rows_selected <- input$dgenes_rows_all
    #   
    # })
    
    
    
    
    
    # observeEvent(input$start_analysis, {
    #   proxy <- dataTableProxy('dgenes', session = session)
    #   selectRows(proxy, selected = NULL)
    #   # rows_selected2 <- input$dgenes_rows_all
    #   rows_selected2 <- NULL
    #   test311113 <<- rows_selected2
    #   rows_selected <- rows_selected2
    # })
    
    
    
    
    tmp_df <- data_selected() %>% 
      select(-start, -end, -chrom) %>%
      filter(source == 'CNV')
    
    
    test777777 <<- rows_selected
    
    if (length(rows_selected) == nrow(tmp_df) | is.null(rows_selected)) {
      tablerInfoCard(
        width = 12,
        value = paste0(nrow(tmp_df), " genes"),
        status = "success",
        icon = "database",
        description = 'Genes found in CNV'
      )
    } else {
      tablerInfoCard(
        width = 12,
        # value = paste0(length(rows_selected), '/', nrow(tmp_df) , " genes"),
        value = paste0(nrow(tmp_df), " genes"),
        
        status = "success",
        icon = "database",
        description = 'Genes found in CNV'
      )
    }
    
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
    
    req(input$input_geno_karyo != 'Multiple coordinates')
    
    
    
    if (input$input_geno_karyo == 'Genomic coordinates') {
      
      
      tmp_cyto <- chromPlot::hg_cytoBandIdeo %>%
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
      description =  'Cytoband(s) selected'
      
    )
    
    
  })
  
  
  output$n_genes_enh_added <- renderUI({
    
    req(nrow(data_selected_enhancers()) > 0)

    tmp_df <- data_selected_enhancers()
    
        tablerInfoCard(
      width = 12,
      value = paste0('+', nrow(tmp_df), " genes"),
      status = "warning",
      icon = "database",
      # description =  name_region
      description = 'Target-genes enhancers'
    )

    
  })
  
  output$n_genes_mirna_added <- renderUI({
    
    req(nrow(data_selected_mirnas()) > 0)
    
    tmp_df <- data_selected_mirnas()
    
    tablerInfoCard(
      width = 12,
      value = paste0('+', nrow(tmp_df), " genes"),
      status = "warning",
      icon = "database",
      # description =  name_region
      description = 'Target-genes miRNAs'
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
    
    req(input$start_analysis > 0)

    data_tmp <- mirtarbase %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap)
    
    
    
  })
  
  
  
  output$df_mirna <- renderDataTable({
    
    
    
    tmp_df <- mirna_raw() %>%
      mutate(references = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', references),
                                 "' target='_blank'>", references,"</a>")) %>%
      mutate(name = paste0("<a href='", paste0('http://www.mirbase.org/textsearch.shtml?q=', name),
                           "' target='_blank'>", name,"</a>")) %>%
      mutate(id = paste0("<a href='", paste0('http://mirtarbase.cuhk.edu.cn/php/detail.php?mirtid=', id),
                         "' target='_blank'>", id,"</a>")) %>%
      select(-contains('source'))
    

    datatable(tmp_df, rownames = FALSE, escape = FALSE,
              colnames = c('ID', 'Name', 'Chrom', 'Start', 'End',
                           'Target-gene', 'Validation experiment', 'Reference'))
    
    
    
  })
  
  tf_raw <- reactive({
    req(input$start_analysis > 0)

    test99993 <<- coord_user()
    
    data_tmp <- trrust %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      select(-target_chrom, -target_start, -target_end) %>%
      mutate(reference = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', reference),"' target='_blank'>", reference,"</a>"))
    
    data_tmp
    
  })
  
  
  output$tf_df <- renderDataTable({
    
    tmp_tbl <<- tf_raw() %>% select(-contains('source'))
    
    datatable(tmp_tbl, rownames = FALSE, escape = FALSE,  
              colnames = c('TF', 'Chrom', 'Start', 'End', 'Target-gene', 'Mechanism', 'Reference'))
    
    
    
  })
  
  
  lncrna_raw <- reactive({
    
    req(input$start_analysis > 0)

    
    data_tmp <- lncrna_coord %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      select(id) %>%
      distinct() %>%
      pull(id)
    
    
    
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
    
    req(input$start_analysis > 0)
    
    data_tmp <- df_enhancers %>% 
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap)
    
    
    # data_tmp <- data_tmp %>% 
    #   left_join(hgcn_genes %>% select(gene, chrom, start, end) %>% 
    #               rename(chrom_gene = chrom), by = 'gene') %>%
    #   rowwise() %>%
    #   mutate(inside_cnv = c(start_gene, end_gene) %overlaps% c(start_coordinates, end_coordinates)) %>%
    #   mutate(inside_cnv = if_else(chrom_gene == chrom_coordinates, inside_cnv, FALSE)) %>%
    #   select(-chrom_gene, -start_gene, -end_gene) %>%
    #   distinct() %>%
    #   filter(!is.na(inside_cnv))
    
    # test1946 <<- data_tmp

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
      title = "Enhancers",
      # trend = -10,
      width = 12
    )
    
    
    # tablerInfoCard(
    #   value =  length(data_tmp),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of Enhancers disrupted",
    #   width = 12
    # )
    
  })
  
  
  # redundancy_enhancers <- reactive({
  #   
  #   data_tmp <- prev_enhancer() %>% select(gene) %>% distinct() %>% pull(gene)
  #   df_output <- df_enhancers %>% filter(gene %in% data_tmp) %>% count(gene)
  #   test25 <<- df_output
  # })
  
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
  
  output$p100_enhancer <- renderPlot({
    
    validate(
      need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
    )
    
    score_filtered <- prev_enhancer() %>% slice(input$df_enhancer_rows_selected) %>% select(phast100) %>% pull()
    
    # prev_enhancer() %>% ggplot(aes(phast100)) +
    #   geom_histogram() +
    #   theme_fancy() +
    #   geom_vline(xintercept = score_filtered, color = 'red')
    plot_p100 +
      theme_fancy() +
      geom_vline(xintercept = score_filtered, color = 'red')
    
    
    
  })
  
  output$p46_enhancer <- renderPlot({
    
    validate(
      need(input$df_enhancer_rows_selected != '', "Please, select an enhancer in the datatable.")
    )
    
    score_filtered <- prev_enhancer() %>% slice(input$df_enhancer_rows_selected) %>% select(phast100) %>% pull()
    
    
    plot_p46pla +
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

    df_tmp <- prev_enhancer()  %>% select(-score_enh, -score) %>%
      select(id,chrom, start, end, gene, everything())
      # filter(!is.na(inside_cnv))
    datatable(df_tmp, rownames = FALSE, filter = 'top', selection = 'single',
              colnames = c('ID Enhancer', 'Chrom',
                           'Start',  'End','Target-gene', 'Phast100way', 'Phast46way Placental', 'Phast46way Primates','o/e gnomad'),
              options = list(
                pageLength = 5, 
                # autoWidth = TRUE,
                # style = 'bootstrap',
                list(searchHighlight = TRUE)
                # selection = 'single'
                # columnDefs = list(list(className = 'dt-center', targets = '_all'))
              ))
    
    
  })
  
  number_tads <- reactive({

    req(input$start_analysis > 0)

    n_tads <- check_tads(coord_user(), tad )
    n_tads <- nrow(n_tads)
    
    if (is.null(n_tads)) {
      n_tads <- 0
    } else {
      n_tads <- as.double(n_tads)
    }
    n_tads

  })
  
  output$n_tads <- renderUI({
    
    
    tablerStatCard(
      value =  number_tads(),
      title = "TADs",
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  number_tads(),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Number of TADs disrupted",
    #   width = 12
    # )
    
    
  })
  
  
  output$n_filtered_enhancers <- renderUI({
    
    
    if (is.null(input$df_enhancer_rows_all)) {
      n_enhancers <- prev_enhancer()
    } else {
      n_enhancers <- prev_enhancer()[input$df_enhancer_rows_all,]
    }
    
    n_total_enhancers <- prev_enhancer() %>% select(gene) %>% distinct() %>% nrow()
    n_enhancers <- n_enhancers %>% select(gene) %>% distinct() %>% nrow()
    
    
    tablerInfoCard(
      width = 12,
      value =  paste0(n_enhancers, '/', n_total_enhancers, ' target-genes'),
      status = "warning",
      icon = "database"
      # description =  ''
      
    )
    
  })
  
  output$switch_tads <- renderUI({
    
    
    validate(
      need(number_tads() > 0, "Need a region with at least one enhancer.")
    )
    
    
    switchInput(
      inputId = "tads_on_off",
      label = "Add genes to analysis?",
      inline = TRUE,
      width = 'auto',
      # status = "warning",
      value = FALSE
      # right = TRUE
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
  
  output$switch_mirnas <- renderUI({
    
    switchInput(
      inputId = "mirnas_on_off",
      label = "Add target genes to analysis?",
      inline = TRUE,
      width = 'auto',
      value = FALSE
    )
    
  })
  
  
  
  prev_tads <- reactive({
    
    req(input$start_analysis > 0)
    

    tmp_df <- data_selected_prev()
    

    
    n_tads <- check_tads(coord_user(), tad )
    
    validate(
      need(!is.null(nrow(n_tads)), "0 TADs found.")
    )
    
    n_tads <- n_tads %>%
      mutate(n_genes = NA) %>%
      mutate(n_genes_not_cnv = NA) %>%
      mutate(genes_not_cnv = NA)
    
    
    for (i in 1:nrow(n_tads)){
      
      tmp_start <- n_tads$start[i]
      tmp_end <- n_tads$end[i]
      
      tmp_genes <- hgcn_genes %>%
        bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
        select(-startdelete, -enddelete, -.overlap) %>%
        pull(gene)
      
      genes_not_cnv <- tmp_genes[!tmp_genes %in% (tmp_df %>% pull(gene))]
      
      n_tads$n_genes[i] <- length(tmp_genes)
      n_tads$n_genes_not_cnv[i] <- length(genes_not_cnv)
      n_tads$genes_not_cnv[i] <- paste(genes_not_cnv, collapse = ', ')
      
    }
    n_tads
  })
  
  
  output$df_tads <- renderDataTable({
    
    tmp_df <- prev_tads() %>% select(-n_genes_not_cnv, -genes_not_cnv)
    
    
    datatable(tmp_df, 
              rownames= FALSE,
              colnames = c('ID', 'Chromosome', 'Start', 'End', 'Total nº of genes')
    )
    
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
    
    # tablerInfoCard(
    #   value = paste(n_dev_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Developmental disorder genes",
    #   width = 12
    # )
    
  })
  
  output$n_clinvar <- renderUI({
    
    
    tablerStatCard(
      value =  nrow(running_clinvar()),
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
    
    # tablerInfoCard(
    #   value =  paste(n_disease_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description = "Disease genes",
    #   width = 12
    # )
    
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
  
  output$hpo_unique_genes_panel <- renderUI({
    
    hpo_yes <-  data_selected() %>% select(gene) %>% pull()
    
    n_genes <- hpo_genes %>% 
      filter(gene %in% hpo_yes) %>% 
      select(gene) %>% 
      distinct() %>%
      nrow()
    
    tablerInfoCard(
      width = 12,
      value =  paste(n_genes, 'Genes'),
      status = "primary",
      icon = "database",
      description =  'Genes with HPO terms'
      
    )
    
    
    
  })
  
  
  hpo_filter <- reactive({
    
    hpo_yes <-  data_selected() %>% select(gene) %>% pull()
    
    hpo_genes_filter <- hpo_genes %>% 
      filter(gene %in% hpo_yes)

  })
  
  
  # output$comparison_patho_pheno <- renderPlotly({
  #   
  # 
  #   p <- data_selected() %>% inner_join(test_hpo(), by = 'gene') %>%
  #     select(gene, pLI, similarity_score) %>%
  #     ggplot(aes(pLI, similarity_score)) +
  #     geom_point(shape = 21, aes(text = gene), fill = 'steelblue', color = 'black', show.legend = FALSE) +
  #     xlim(c(0,100)) +
  #     ylab('Similarity score') +
  #     xlab('pLI score (0-100)') +
  #     theme_fancy()
  #   
  #   ggplotly(p)
  #   
  # 
  #   
  # })
  # 
  
  # output$hpo_filter_genes <- renderDataTable({
  #   
  #   
  #   # validate(
  #   #   need(length(input$chosen_hp) > 0, "Please, select at least one HPO term.")
  #   # )
  #   
  #   # validate(
  #   #   need(input$run_pheno_analysis, "Click on run analysis.")
  #   #   
  #   # )
  #   
  #   hpo_yes <-  data_selected() %>% select(gene) %>% pull()
  #   
  #   n_genes <- hpo_genes %>% 
  #     filter(gene %in% hpo_yes) %>% 
  #     select(gene) %>% 
  #     distinct() %>%
  #     nrow()
  #   
  #   validate(
  #     need(n_genes != 0, "0 genes associated with HPO terms.")
  #     
  #   )
  #   
  #   tmp_df <- hpo_filter()   %>% 
  #     count(gene) %>% 
  #     left_join(test_hpo(), by = 'gene') %>%
  #     left_join((data_selected() %>% select(source, gene)), by = 'gene') %>%
  #     select(source, everything()) %>%
  #     mutate(similarity_score = replace_na(similarity_score, '-'))
  #   
  #   datatable(tmp_df, selection = 'single', rownames = FALSE,
  #             filter = 'top', colnames = c('Source','Gene', 'Nº HPO terms', 'Similarity score')
  #             
  #   )
  #   
  # })
  
  
  running_sim_score <- reactive({
    
    req(input$start_analysis > 0)
    
    req(input$input_inheritance != '')
    
    
    tmp_tbl <- test2019 %>% 
      select(source, gene) %>%
      left_join(hpo_genes %>% select(identifier, gene, hp), by = 'gene') %>%
      rename(hp_gene = hp) %>%
      na.omit()
    
    if (input$input_inheritance != 'Any') {
    keep_identifiers <- tmp_tbl %>% select(identifier) %>%
      left_join(hpo_omim %>% select(identifier, hp)) %>%
      filter(hp %in% input$input_inheritance) %>%
      pull(identifier)
    
    tmp_tbl <- tmp_tbl %>% filter(identifier %in% keep_identifiers)
    
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
      
      # hpo_patient  <- list('patient' = 'HP:0001285')
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
    
    output$dt_running_sim_score <- renderDT({
      
      
      datatable(running_sim_score() %>% replace_na(list(sim_gene = '-', sim_mim = '-')) %>%
                  mutate(identifier = paste0(identifier, '(', disease_source, ')')) %>%
                  select(-disease_source), 
                rownames = FALSE,
                selection = 'single',
                colnames = c('Source', 'Gene', 'Identifier', 'Disease description', 'Similarity score (gene)',
                             'Similarity score (disease)'))
    })

  
  
  # run_sim_diseases <- reactive({
  #   
  #   tmp <- data_selected() %>% select(gene) %>% pull()
  # 
  #   input_mim_disease <- omim %>%
  #     filter(gene %in% tmp) %>%
  #     select(MIM_pheno_number, gene) %>%
  #     rename(mim_disease = MIM_pheno_number)
  # 
  #   n_hpo_terms <- hpo_omim %>%
  #     count(mim_disease, desc, name = "n_hpo_terms" ) %>%
  #     inner_join(input_mim_disease, by = 'mim_disease')
  #   
  #   validate(
  #     need(nrow(n_hpo_terms) != 0, "0 OMIM diseases associated with HPO terms.")
  #     
  #   )
  #   
  #   # validate(
  #   #   need(length(input$chosen_hp) > 0, "Please, select at least one HPO term.")
  #   # )
  #   
  #   # validate(
  #   #   need(input$run_pheno_analysis, "Click on run analysis.")
  #   #   
  #   # )
  #   
  #   tmp_df <- n_hpo_terms %>%
  #     left_join(calculate_sim_disease(), by = 'mim_disease') %>%
  #     mutate(mim_disease = paste0("<a href='", paste0('https://www.omim.org/entry/', mim_disease),"' target='_blank'>", mim_disease,"</a>")) %>%
  #     select(-n_hpo_terms) %>%
  #     mutate(similarity_score = replace_na(similarity_score, '-'))
  #   
  # })
  
  
  output$n_diseases <- renderUI({
    
    
    number_diseases <- running_sim_score() %>%  nrow()

    
    tablerInfoCard(
      width = 12,
      value =  paste(number_diseases, 'diseases'),
      status = "primary",
      icon = "database",
      description =  ''
    )
    
  })
  
  output$hpo_filter_diseases <- renderDataTable({
    
    
    datatable(run_sim_diseases(), 
              rownames = FALSE, 
              escape = FALSE, 
              selection = 'single', 
              colnames = c('MIM term', 'Disease Name','Gene', 'Similarity score'))
    
    
  })
  
  
  

  
  # output$hpo_assoc_genes <- renderDataTable({
  #   
  #   # hpo_filter_genes
  #   
  #   validate(
  #     need(input$hpo_filter_genes_rows_selected != '', "Please, select a gene.")
  #   )
  #   
  #   
  #   
  #   tmp_df <- hpo_filter() %>%
  #     count(gene) %>%
  #     slice(input$hpo_filter_genes_rows_selected) %>%
  #     select(gene) %>%
  #     pull(gene)
  #   
  #   
  #   tmp_df2 <- hpo_filter()
  #   tmp_df2 <- tmp_df2 %>% filter(gene %in% tmp_df ) %>% select(-entrez_id)
  #   tmp_df2 <- tmp_df2 %>% 
  #     mutate(hp = paste0("<a href='", paste0('https://hpo.jax.org/app/browse/term/', hp),"' target='_blank'>", hp,"</a>")) %>%
  #     select(gene, hp, term)
  #   test232134 <<- tmp_df2
  #   datatable(tmp_df2, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'HPO term', 'Description'))
  #   
  # })
  
  
  
  # run_select_disease_hpo <- reactive({
  #   
  #   
  # 
  #   
  # })
  # 
  output$hpo_assoc_diseases <- renderDataTable({
    
    validate(
      need(input$dt_running_sim_score_rows_selected != '', "Please, select a row.")
    )
    
    test301 <<- running_sim_score()
    
    
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
                                         paste0('https://www.orpha.net/consor/cgi-bin/Disease_Search.php?lng=EN&data_id=8648&Disease_Disease_Search_diseaseGroup=', identifier),"' target='_blank'>", identifier,"</a>"))
      ) %>%
      mutate(identifier = paste0(identifier, ' (', disease_source, ')')) %>%
      select(-disease_source)
                                    
    
    
    datatable(tmp_df, escape = FALSE, rownames = FALSE, colnames = c('Disease identifier', 'HPO term', 'Description'))
    
  })
  
  output$hpo_assoc_genes <- renderDataTable({
    
    validate(
      need(input$dt_running_sim_score_rows_selected != '', "Please, select a row.")
    )
    
    tmp_df <- running_sim_score() %>%
      slice(input$dt_running_sim_score_rows_selected) %>%
      pull(gene)
    
    
    test91212411 <<- tmp_df
    
    tmp_df2 <- hpo_filter() %>% 
      filter(gene %in% tmp_df)
      # left_join(hpo_genes %>% select(term, hp) %>% distinct(), by = c('hp' = 'term')) 
    
    test9999999999999999 <<- tmp_df2
    
    tmp_df2 <- tmp_df2 %>% 
      mutate(hp = paste0("<a href='", paste0('https://hpo.jax.org/app/browse/term/', hp),"' target='_blank'>", hp,"</a>")) %>%
      select(gene, hp, desc)
    
    datatable(tmp_df2, escape = FALSE, rownames = FALSE, colnames = c('Gene', 'HPO term', 'Description'))
    
  })
  
  
  # output$hpo_filter_cnvs <- renderDataTable({
  #   
  #   
  #   test913 <<- data_selected()
  #   test914 <<- test_hpo()
  #   
  #   tmp_df <- hpo_filter() 
  #   
  #   tmp_df <- tmp_df  %>% count(gene) %>% left_join(test_hpo(), by = 'gene') 
  #   datatable(tmp_df, filter = 'top', colnames = c('Gene', 'Nº HPO terms', 'Similarity score'), option = list(
  #     selection = 'single'
  #   ))
  #   
  # })
  
  
  # output$gene_filter_hp <- renderUI({
  #   
  #   test1452 <<- test_hpo()
  #   
  #   max_value <- test_hpo() %>% select(similarity_score) %>%
  #     arrange(desc(similarity_score)) %>%
  #     slice(1) %>%
  #     pull()
  #   
  #   test1451 <<- max_value
  #   
  #   sliderTextInput(
  #     inputId = "input_gene_filter_hp",
  #     label = "Choose a value:", 
  #     choices = round(seq(from = 0, to = 20, by = 1)),
  #     grid = TRUE
  #   )
  #   
  #   
  #   
  #   
  #   
  # })
  
  output$plot_anatomy <- renderPlot({
    
    
    validate(
        need(length(input$chosen_hp) != 0, 'Please, select at least one HPO term.'),
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
    
    
    test100 <<- hpo_from_gene
    test200 <<- hpo_from_patient
    test300 <<- hpo_from_disease
    
    
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
    
    
    test91312 <<- plot_df
    
    plot_df %>%
      mutate(value = as.factor(value)) %>%
      mutate(value = fct_relevel(value, 'Ear', after = 0 )) %>%
      mutate(value = fct_relevel(value, 'Breast', after = 1 )) %>%
      ggplot(aes(x = value, y = valuae)) +
      geom_col(aes(fill = class), position = 'dodge', color = 'black') +
      theme_fancy() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            legend.position = 'top',
            legend.text=element_text(size= 13)) +
      labs( y = 'Number of HPO terms', 
            x = 'Anatomical entities', 
            fill = NULL)
    
    
    
  })
  

  output$plot_similarity_genes <- renderPlot({
    
    
    validate(
      need(length(input$chosen_hp) > 0, "Please, select at least one HPO term.")
    )
    
    
    tmp_df <- running_sim_score()
    

    if (input$select_sim_gene_disease == 'genes') {

      tmp_df  %>%
        select(gene, sim_gene) %>%
        arrange(desc(sim_gene)) %>%
        distinct() %>%
        ggplot(aes(reorder(gene, sim_gene), sim_gene)) +
        geom_col(aes(fill = sim_gene), color = 'black', show.legend = FALSE) +
        ylab('Phenotypic similarity score') +
        xlab('Genes') +
        scale_fill_viridis_c() +
        theme_fancy() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16))

    } else {

      tmp_df  %>%
        select(identifier, sim_mim) %>%
        arrange(desc(sim_mim)) %>%
        # filter(similarity_score >= filter_higher_than) %>%
        ggplot(aes(reorder(identifier, sim_mim), sim_mim)) +
        geom_col(aes(fill = sim_mim), color = 'black', show.legend = FALSE) +
        # coord_flip() +
        ylab('Phenotypic similarity score') +
        xlab('Disease identifiers') +
        scale_fill_viridis_c() +
        theme_fancy() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16))

    }


  })
  
  
  output$n_hi<- renderUI({
    
    
    n_hi_yes <- data_selected() %>% filter(hi <= 10) %>% nrow()
    n_total <- nrow(data_selected())  
    
    tablerStatCard(
      value =  paste(n_hi_yes, n_total, sep = '/'),
      title =  paste("Likely to exhibit haploinsufficiency"),
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  paste(n_hi_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description =paste("Likely to exhibit haploinsufficiency", '(hi <= 10)'),
    #   width = 12
    # )
  })
  
  output$n_pli <- renderUI({
    
    test15 <- data_selected() %>% select(entrez_id) %>% pull()
    
    
    n_pli_yes <- data_selected() %>% filter(pLI >= 95) %>% nrow()
    n_total <- nrow(data_selected())  
    
    tablerStatCard(
      value =  paste(n_pli_yes, n_total, sep = '/'),
      title =  paste("Intolerant to LoF mutations", '(pLI >= 0.9)'),
      # trend = -10,
      width = 12
    )
    
    # tablerInfoCard(
    #   value =  paste(n_pli_yes, n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description =  paste("Intolerant to LoF mutations", '(pLI >= 0.9)'),
    #   width = 12
    # )
    
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
    
    
    # tablerInfoCard(
    #   value =  paste(length(gene_id), n_total, sep = '/'),
    #   status = "primary",
    #   icon = 'book',
    #   description = "OMIM genes",
    #   width = 12
    # )
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
    
    
    test_tmp <- data_selected() %>% select(gene) %>% pull()
    mgi_tmp <- mgi %>% filter(gene %in% test_tmp)
    mgi_tmp <- mgi_tmp %>% separate(pheno, into = LETTERS[1:230], sep = ' ') %>%
      gather('delete', 'mpo_id', -gene, -entrez_id, -gene_mouse, -mgi) %>%
      filter(mpo_id != '') %>%
      select(-delete)
    
    test2141241241 <<- mgi_tmp
    
    vector_mpo <- mgi_tmp %>% 
      select(mpo_id) %>% 
      distinct() %>%
      left_join(mpo_dbs, by = c('mpo_id' = 'term'))
    
    mgi_tmp <- mgi_tmp %>% left_join(vector_mpo, by = 'mpo_id')
    
    mgi_tmp
    
    
    
  })
  
  
  output$agg_model <- renderPlot({
    
    
    tmp_df <- model_genes_phenotype() %>% 
      count(description) %>% arrange(desc(n))
    
    if (nrow(tmp_df) > 10) {
      
      tmp_df <- tmp_df %>% slice(1:10)
    }
    
    tmp_df %>%
      ggplot(aes(reorder(description, n), n)) + 
      geom_col(aes(fill = n), color = 'black') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      coord_flip() +
      scale_fill_viridis_c() +
      ylab('Number of genes') +
      xlab('Mouse phenotype terms') + 
      theme_fancy() +
      theme(axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16)) +
      labs(fill = NULL) +
      scale_y_continuous(breaks = scales::pretty_breaks())
    
    # model_genes_phenotype() %>% count(description) %>% arrange(n) %>%
    #   e_charts() %>% 
    #   e_treemap(description, description, n) %>%
    #   e_tooltip(trigger = "axis") %>%
    #   e_title("")
    
    
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
    test231312321 <<- go_analysis
    
    
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
    
    test01441414 <<- filtered_genes
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
    
    test9991 <<- running_path()
    df <- test9991  %>%
      # filter(Count != 0) %>%
      arrange(desc(Count)) %>%
      separate(geneID, sep = '/', into = as.character(1:1000)) %>%
      gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
      select(-delete) %>%
      # na.omit() %>%
      distinct()
    
    datatable(df, escape = FALSE)
    
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
      need(length(filtered_genes) != 0, "0 enriched diseases found.")
    )
    
    enrich_dgn <- enrichDO(gene  = filtered_genes,
                           universe      = hgcn_genes %>% select(entrez_id) %>% pull() %>% as.character(),
                           # pAdjustMethod = "BH",
                           pvalueCutoff  = as.numeric(input$pvalue_do),
                           readable = TRUE)
    
    test444 <<- enrich_dgn %>% as_tibble() %>% filter(n > 1)
    
    validate(
      need(nrow(enrich_dgn) != 0, "0 enriched terms found.")
    )
    
    enrich_dgn
    
  })
  
  
  output$df_do <- renderDataTable({
    
    test883 <<- running_do()  
    
    df <- running_do() %>%
      as_tibble() %>%
      filter(Count != 0) %>%
      arrange(desc(Count)) %>%
      separate(geneID, sep = '/', into = as.character(1:1000)) %>%
      gather('delete', 'gene', -ID, -Description, -Count, -GeneRatio, -pvalue, -p.adjust, -qvalue, -BgRatio) %>%
      select(-delete, -qvalue) %>%
      na.omit() %>%
      distinct() %>%
      mutate(pvalue = round(pvalue, 3)) %>%
      mutate(p.adjust = round(p.adjust, 3)) %>%
      mutate(qvalue = round(qvalue, 3)) %>%
      mutate(ID = paste0("<a href='", paste0('https://www.ebi.ac.uk/ols/ontologies/doid/terms?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F', ID),"' target='_blank'>", ID,"</a>"))
    
    
    datatable(df, escape = FALSE)
    
  })
  
  output$func_do  <- renderPlot({
    
    
    cnetplot(running_do(), foldChange= hgcn_genes$entrez_id, readable = TRUE)
  })
  
  # output$table_all_genes  <-  renderReactable({
  #   
  #  reactable(hgcn_genes,
  #            filterable = TRUE,
  #            selection = "single", 
  #            selectionId = "selected",
  #            
  #            columns = list(chrom = colDef(name = "Chromosome"),
  #                           start_position = colDef(name = "Start"),
  #                           end_position = colDef(name = "End"),
  #                           location = colDef(name = "Location"),
  #                           gene = colDef(name = "Gene"),
  #                           haplo = colDef(name = "Haploinsufficiency",
  #                                          cell = function(value) {
  #                                            
  #                                            if (value == 0) "\u2718" else "\u2713"
  #                                          }),
  #                           triplo = colDef(name = "Triplosensitivity",
  #                                           cell = function(value) {
  #                                             
  #                                             if (value == 0) "\u2718" else "\u2713"
  #                                           }),
  #                           dev = colDef(name = "Developmental disorder gene",
  #                                        cell = function(value) {
  #                                          
  #                                          if (value == 0) "\u2718" else "\u2713"
  #                                        }),
  #                           fda = colDef(name = "FDA-approved drug targets",
  #                                        cell = function(value) {
  #                                          
  #                                          if (value == 0) "\u2718" else "\u2713"
  #                                        }),
  #                           clinvar = colDef(name = "ClinVar genes",
  #                                            cell = function(value) {
  #                                              
  #                                              if (value == 0) "\u2718" else "\u2713"
  #                                            }),
  #                           ccr = colDef(name = "CCR"),
  #                           ncrvis = colDef(name = "Ncrvis"),
  #                           ncgerp = colDef(name = "Ncgerp"),
  #                           hi = colDef(name = "HI index", 
  #                                       format = colFormat(suffix = " %", digits = 1))
  #                           # cell = function(value) {
  #                           #   if (is.na(value)) {
  #                           #     classes <- "tag num-low"
  #                           #   } else if (value >= 0.9) {
  #                           #     classes <- 'tag num-high'
  #                           #   } else  {
  #                           #     classes <- "tag num-high"
  #                           #   }
  #                           #   value <- format(value, nsmall = 1)
  #                           #   span(class = classes, value)
  #                           # })
  #            ))
  #   
  # })
  
  
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

    req(input$input_geno_karyo == 'Multiple coordinates')
    
      validate(
        need(!is.null(input$file_cnv), "Please upload a file.")
      )

    file1 <- input$file_cnv
    test0020 <<- file1
    data1 <- read_tsv(file1$datapath, col_names = c('chrom', 'start', 'end'),
                      col_types = list(chrom = col_character(),
                                       start = col_integer(),
                                       end = col_integer()))
    data1 
  })
  
  output$cnv_file <- renderDT({
    
    tmp_tbl <- reading_cnv_file() %>% mutate(length_cnv = end - start + 1)
    
    datatable(tmp_tbl, colnames = c('Chrom', 'Start', 'End', 'Lenght'))
    
  })
  
  cnv_file_to_analyze <- reactive({
    
    tmp_tbl <- reading_cnv_file()
    

    if (input$select_n_cnvs == 'yes') {
      
      tmp_tbl
      
    } else {
      
      tmp_tbl[input$cnv_file_rows_all,]
    }
    
    
    
    
    
    
  })
  
  df_overlap_cnvs_running <- reactive({
    
    req(input$start_analysis > 0)

    tmp_df <-  get_perc_overlap(check_cnv_df(),
                                coord_user())
    
    
    tmp_df
  })
  
  # intersection_running <- reactive({
  # 
  #   start_coordinates <- as.numeric(coord_user()[1])
  #   end_coordinates <- as.numeric(coord_user()[2])
  #   chrom_coordinates <- coord_user()[3]
  #   
  #   test00000 <<- df_overlap_cnvs_running() 
  # 
  #   df_pathogenic <- df_overlap_cnvs_running() %>% filter(source == 'decipher') 
  #   
  #   test93131 <<- df_pathogenic
  #   
  #   validate(
  #     need(nrow(df_pathogenic) != 0, "No pathogenic CNVs found")
  #   )
  # 
  #   cnv_input <- GRanges(
  #     seqnames= chrom_coordinates,
  #     ranges= IRanges(start= start_coordinates, end = end_coordinates)
  #   )
  # 
  #   cnvs_pathogenic <- GRanges(
  #     seqnames= df_pathogenic$chrom,
  #     ranges= IRanges(start= df_pathogenic$start_position, end = df_pathogenic$end_position)
  #   )
  # 
  #   gr_intersect <- intersect(cnv_input, cnvs_pathogenic) %>% as_tibble()
  #   
  #   gr_intersect <- gr_intersect %>%
  #     mutate(id = as.factor(seq(1:n())))
  # 
  #   test711 <<- gr_intersect
  # 
  #   gr_intersect
  # 
  # 
  # })
  
  output$df_intersection <- renderDataTable({
    
    
    df_tmp <- intersection_running() %>%
      select(id, start, end)
    
    datatable(df_tmp, rownames = FALSE, selection = 'single',
              options = list(dom = 't'))
    
  })
  
  # output$plot_intersection <- renderPlot({
  #   
  #   validate(
  #     need(FALSE, "TBA")
  #   )
  #   
  # 
  #   start_coordinates <- as.numeric(coord_user()[1])
  #   end_coordinates <- as.numeric(coord_user()[2])
  #   
  # 
  #   
  #   df_tmp <-intersection_running()
  # 
  #   df_genes <- data_selected() %>% mutate(keep = NA) %>%
  #     rowwise() %>%
  #     mutate(keep = c(start_position, end_position) %overlaps% c(df_tmp$start, df_tmp$end)) %>%
  #     ungroup() %>%
  #     filter(keep == TRUE) %>%
  #     select(-keep) %>%
  #     ungroup() %>%
  #     mutate(pos = round((end_position + start_position)/2), 0) %>%
  #     mutate(id = rep('gene(s)', n())) %>%
  #     select(id, pos)
  # 
  #   df_total <- intersection_running() %>%
  #     select(id, start, end) %>%
  #     gather('rm', 'pos', -id) %>%
  #     select(-rm) %>%
  #     rbind(df_genes) %>%
  #     mutate(id_color = as.factor(if_else(id == 'gene(s)', 'steelblue', 'red')))
  #   
  #     test1453 <<- df_total
  # 
  #   ggplot() +
  #     geom_point(data = df_total, aes(pos, id, fill = id_color), color = 'black', shape = 21) +
  #     geom_path(data =df_total %>% filter(id != 'gene(s)'), aes(pos, id, group = id, color = id_color)) +
  #     coord_cartesian(xlim = c(start_coordinates, end_coordinates )) +
  #     ylab('Id intersections') +
  #     xlab(paste0('chr', coord_user()[3], ':', coord_user()[1], '-', coord_user()[2])) +
  #     theme_fancy()
  # 
  # })
  
  
  df_overlap_cnvs_running_download_patho <- reactive({
    
    
    tmp_df <-  df_overlap_cnvs_running() %>% filter(source == 'decipher')
    
    tmp_df
  })
  
  df_overlap_cnvs_running_download_nonpatho <- reactive({
    
    
    tmp_df <-  df_overlap_cnvs_running() %>% filter(source != 'decipher')
    
    tmp_df
  })
  
  # output$plot_cnv_bar <- renderPlot({
  #   
  #   
  #   start_pos <- as.numeric(coord_user()[1])
  #   end_pos <- as.numeric(coord_user()[2])  
  #   
  #   
  #   
  #   df_nonpathogenic <- df_overlap_cnvs_running() %>% filter(source != 'decipher')
  #   df_pathogenic <- df_overlap_cnvs_running() %>% filter(source == 'decipher')
  #   
  #   test8231 <<- start_pos
  #   test8999 <<- end_pos
  #   test88881 <<- df_nonpathogenic
  #   test88882 <<- df_pathogenic
  #   
  #   # start_pos <- test8231
  #   # end_pos <- test8999
  #   # df_nonpathogenic <- test88881
  #   # df_pathogenic <- test88882
  #   
  #   
  #   interval_values <- seq(from = start_pos, to = end_pos, length.out = 200)
  #   
  #   df_interval <- matrix(interval_values, ncol = 2, byrow = TRUE)
  #   colnames(df_interval) <- c('start', 'end')
  #   df_interval <- as_tibble(df_interval)
  #   
  #   query <- IRanges(df_interval$start, df_interval$end)
  #   
  #   # Evaluate pathogenic CNVs
  #   subject_patho <- IRanges(df_pathogenic$start_position, df_pathogenic$end_position)
  #   hits_intervals_patho <- countOverlaps(query, subject_patho)
  #   
  #   # Evaluate pathogenic CNVs
  #   subject_nonpatho<- IRanges(df_nonpathogenic$start_position, df_nonpathogenic$end_position)
  #   hits_intervals_nonpatho <- countOverlaps(query, subject_nonpatho)
  #   
  #   
  #   df_interval <- df_interval %>% mutate(n_patho = hits_intervals_patho, 
  #                                         n_nonpatho = hits_intervals_nonpatho) %>%
  #     mutate(id = row_number()) %>%
  #     gather('category', 'n_overlap', -start, -end, -id)
  #   
  #   
  #   df_interval %>%
  #     mutate(category = if_else(category == 'n_patho', 'Pathogenic CNVs', 'Non-pathogenic CNVs')) %>%
  #     ggplot(aes(id, n_overlap)) +
  #     geom_col(aes(fill = category), color = 'black', show.legend = FALSE) + 
  #     #geom_point() +
  #     # geom_smooth(aes(color = category, group = category))
  #     #geom_line(aes(color = category, group = category)) + 
  #     theme_minimal() + 
  #     facet_wrap(~category, nrow = 2)
  #   
  # })
  
  
  
  
  output$df_overlap_cnvs <- renderDT({
    

    
    
    
    tmp_df <- df_overlap_cnvs_running() %>% filter(source == 'decipher') %>% 
      filter(pathogenicity %in% c('Pathogenic', 'Likely pathogenic')) %>%
      select(id, chrom, start, end, pathogenicity, genotype, variant_class, phenotypes, length_cnv, p_overlap) %>%
      arrange(desc(p_overlap))
    
    
    validate(
      need(nrow(tmp_df) != 0, "No pathogenic CNVs found.")
    )
    
    datatable(tmp_df, escape = FALSE,
              colnames = c('ID', 'Chrom', 'Start', 'End', 'Pathogenicity', 'Genotype', 'Class', 'Phenotype',
                           'CNV size', 'Overlap (%)'),
              selection = 'single',
              options = list(
                columnDefs = list(list(className = 'dt-center',  targets = c(0:6,8,9)))),  rownames= FALSE)
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
  
  
  output$ui_select_del_dup <- renderUI({
    
    # vector_n_dbs <- df_overlap_cnvs_running()  %>% filter(source != 'decipher') %>%
    #   count(source) %>% 
    #   mutate(new_one = paste0(source, ' (', n, ")")) %>%
    #   pull(new_one)
    # 
    # test142412412 <<- vector_n_dbs
    # 
    # vector_n_dbs <- vector_n_dbs %>% 
    #   str_replace('decipher_control', 'DECIPHER Control') %>% 
    #   str_replace('dgv', 'DGV') %>%
    #   str_replace('gnomad_v2.1', 'gnomAD v.2.1')
    
    vector_n_dbs <-    split(c('deletions', 'duplications'), c('Deletions', 'Duplications'))
    
    prettyRadioButtons(
      inputId = "select_del_dup",
      label = tags$b("Select one option:"),   
      choices =  vector_n_dbs,
      inline = TRUE, 
      status = "primary",
      fill = TRUE
    )
  })
  
  
  
  
  running_dgv <- reactive({
    
    req(input$start_analysis > 0)

    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'dgv') %>% pull(id)
    
    tmp_df <- dgv_df_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      select(id, chrom, start, end, everything())
    
    
    tmp_df
  })
  
  running_decipher_c <- reactive({
    
    req(input$start_analysis > 0)
    
    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'decipher_control') %>% pull(id)
    
    tmp_df <- decipher_control_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      arrange(desc(p_overlap)) %>%
      select(-source) %>%
      select(id, chrom, start, end, everything())
    
    
    
    
    tmp_df
  })
  
  running_gnomad <- reactive({
    
    req(input$start_analysis > 0)

    filter_id <-  df_overlap_cnvs_running() %>%
      filter(source == 'gnomad_v2.1') %>% pull(id)
    
    tmp_df <- gnomad_sv_raw %>% 
      filter(id %in% filter_id) %>%
      get_perc_overlap(coord_user()) %>%
      select(id, chrom, start, end, svtype, AF, p_overlap)
    
    
    
    tmp_df
  })
  
  
  
  
  
  
  output$df_overlap_cnvs_nonpatho <- renderDT({
    
    req(input$select_no_patho_cnv)
    
    
    if (input$select_no_patho_cnv == 'dgv') {
      
      datatable(running_dgv(), rownames = FALSE,
                colnames = c('ID', 'Chrom', 'Start', 'End', 'Type', 'Reference', 'PMID',
                             'Method', 'Sample size', 'Observed gains', 'Observed losses', 'Genes',
                             'Overlap (%)'))
      
      
      
    } else if (input$select_no_patho_cnv == 'decipher_control') {
      
      # validate(
      #   need(nrow(running_decipher_c()) != 0, "No non-pathogenic CNVs from DECIPHER Control found.")
      # )
      
      
      datatable(running_decipher_c(), rownames = FALSE,
                colnames = c('ID', 'Chrom', 'Start', 'End', 'Deletion Obs.', 'Deletion Freq.',
                             'Deletion (SE)', 'Duplication Obs.', 'Duplication Freq.', 'Duplication (SE)',
                             'Observations', 'Frequency', 'Standard Error', 'Type','Sample size', 'Study',
                             'Overlap (%)'))
      
      
      
    } else {
      
      datatable(running_gnomad(), 
                rownames = FALSE,
                colnames = c('ID', 'Chrom', 'Start', 'End', 'Type', 'Allele Frequency', 'Overlap (%)'))
      
      
      
    }
    
    
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
  
  
  # 
  # output$df_upload <- renderDataTable({
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
    
    req(input$start_analysis > 0)

    tmp_df <- gwas_variants %>%
      rename(chrom = CHR_ID, start = CHR_POS) %>%
      mutate(end = start) %>%
      bed_intersect(test2020, suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      rename(pos = start) %>%
      select(-end) %>%
      mutate(pubmed_id = str_extract(LINK, '\\d{8}')) %>%
      select(-LINK)
  })
  
  output$df_gwas <- renderDataTable({
    
    tmp_df <- running_gwas() %>%
      mutate(pubmed_id = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pubmed_id),"' target='_blank'>", pubmed_id,"</a>"))
    
    # test0007 <<- tmp_df
    
    datatable(tmp_df, 
              escape = FALSE,
              colnames = c('Chrom', 'Position','Intergenic', 'Disease trait', 'gene', 'Link study'),
              rownames = FALSE)
    
    
  })
  
  
  running_de_novo <- reactive({
    
    req(input$start_analysis > 0)

    tmp_df <- denovo %>%
      rename(start = Position) %>%
      mutate(end = start) %>%
      bed_intersect(test2020, suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      rename(position = start) %>%
      select(-end)
    
    tmp_df
    
  })
  
  running_clinvar <- reactive({
    
    req(input$start_analysis > 0)
    
    

    tmp_df <- clinvar_variants %>%
      mutate(chrom = as.character(chrom)) %>%
      rename(start = pos) %>%
      mutate(end = start) %>%
      bed_intersect(coord_user(), suffix = c('', 'delete')) %>%
      select(-startdelete, -enddelete, -.overlap) %>%
      rename(pos = start) %>%
      select(-end)
    
    tmp_df
    
  })
  
  
  output$ui_select_clinvar_gwas <- renderUI({
    
    req(running_clinvar())
    req(running_gwas())
    
    
    vector_tmp <- c(paste('ClinVar',paste0('(', nrow(running_clinvar()), ')')),
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
  
  output$df_variants <- renderDataTable({
    
    req(input$select_clinvar_gwas)
    
    if (input$select_clinvar_gwas == 'clinvar') {
      
      
      tmp_df <- running_clinvar() %>% 
        mutate(id = paste0("<a href='", paste0('https://www.ncbi.nlm.nih.gov/clinvar/variation/', id),"' target='_blank'>", id,"</a>"))
      
      
      
      datatable(tmp_df, escape = FALSE, colnames = c('Chrom', 'Position','Reference', 'Alternative','Gene','Clinical significance',
                                                     'Disease Identifier', 
                                                     'Disease name', 'Clinvar ID'), 
                rownames = FALSE
      )
      
      
      
    } else {
      
      tmp_df <- running_gwas() %>%
        mutate(pubmed_id = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', pubmed_id),"' target='_blank'>", pubmed_id,"</a>"))
      
      datatable(tmp_df, 
                escape = FALSE,
                colnames = c('Chrom', 'Position','Intergenic', 'Disease trait', 'gene', 'Link study'),
                rownames = FALSE)
      
      
      
      
    }
    
    
    
    
  })
  
  
  
  
  output$df_de_novo <- renderDataTable({
    
    tmp_df <- running_de_novo() %>% 
      
      mutate(PubmedID = paste0("<a href='", paste0('https://pubmed.ncbi.nlm.nih.gov/', PubmedID),"' target='_blank'>", PubmedID,"</a>"))
    
    # tmp_df <- tmp_df %>% 
    
    
    # tmp_df <- running_de_novo() %>% mu
    
    datatable(tmp_df, escape = FALSE, colnames = c('Chrom', 'Position','Gene', 'Phenotype', 'Study name', 'PubmedID', 'Function Class', 
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