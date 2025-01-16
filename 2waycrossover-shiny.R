library(shiny)
library(nlme)
library(dplyr)

ui <- fluidPage(
  titlePanel("两周期实验ABE分析和Bootstrap 工具"),  # 标题
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "上传 CSV 文件", accept = ".csv"),  # 文件上传
      numericInput("n_subjects", "受试者数量 (n_subjects):", value = 10),  # 输入受试者数量
      numericInput("n_bootstrap", "Bootstrap 次数 (n_bootstrap):", value = 100),  # 输入 Bootstrap 次数
      actionButton("run_abe", "运行 ABE 分析"),  # ABE 分析按钮
      actionButton("run_bootstrap", "运行 Bootstrap 分析")  # Bootstrap 分析按钮
    ),
    mainPanel(
      h3("ABE 分析结果"),  # ABE 结果标题
      tableOutput("abe_result"),  # ABE 结果显示
      h3("Bootstrap 分析结果"),  # Bootstrap 结果标题
      tableOutput("bootstrap_result")  # Bootstrap 结果显示
    )
  )
)


server <- function(input, output) {
  # 监听 ABE 分析按钮
  observeEvent(input$run_abe, {
    req(input$file)  # 确保文件已上传
    
    # 读取上传的 CSV 文件
    data <- read.csv(input$file$datapath)
    
    # 调用 ABE 分析函数
    abe_result <- abe_analysis(data)
    
    # 显示 ABE 分析结果
    output$abe_result <- renderTable({
      abe_result
    })
  })
  
  # 监听 Bootstrap 分析按钮
  observeEvent(input$run_bootstrap, {
    req(input$file)  # 确保文件已上传
    
    # 读取上传的 CSV 文件
    data <- read.csv(input$file$datapath)
    
    # 调用 Bootstrap 分析函数
    test_result <- bootstrap_ABE(data, input$n_subjects, input$n_bootstrap)
    
    bootstrap_result<- test_result %>%
      group_by(Parameter) %>%
      summarise(Successcount = sum(success),
                Totalcount = n(),
                success_rate=(Successcount/Totalcount))
    # 显示 Bootstrap 分析结果
    output$bootstrap_result <- renderTable({
      bootstrap_result
    })
  })
}
abe_analysis <- function(data){
  #ABEexample1 <- read.table(file="clipboard", header=TRUE)
  ABEexample1 <- data
  
  #For MacOS:
  #ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)
  
  #Now copy and paste the following into R: (or select "Edit", then "Run all")
  p <- 3 # rounding for Cmax and AUCs
  pci <- 2 # rounding for ratio and 90%CI (use 1 for TPD, 2 for EMA, FDA)
  
  
  factors       <- c("SUBJ", "GRP", "PRD", "TRT")
  ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for lme()
  
  lm.lnauct <- lme(log(AUClast)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  lm.lnaucinf <- lme(log(AUCinf)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  lm.lncmax <- lme(log(Cmax)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  
  #lm.lnauct <- lme(-log(AUClast)~TRT+PRD+GRP+SUBJ, data=ABEexample1,na.action=na.exclude)
  #lm.lnaucinf <- lme(-log(AUCinf)~TRT+PRD+GRP+SUBJ, data=ABEexample1, na.action=na.exclude)
  #lm.lncmax <- lme(-log(Cmax)~TRT+PRD+GRP+SUBJ, data=ABEexample1, na.action=na.exclude)
  
  
  
  #Analysis for AUCt
  summary(lm.lnauct)
  intervals(lm.lnauct,0.90,"fixed")
  VarCorr(lm.lnauct)
  
  #Analysis for AUCinf
  summary(lm.lnaucinf)
  intervals(lm.lnaucinf,0.90,"fixed")
  VarCorr(lm.lnaucinf)
  
  #Analysis for Cmax
  summary(lm.lncmax)
  intervals(lm.lncmax,0.90,"fixed")
  VarCorr(lm.lncmax)
  
  #roll up the answers
  Parameter <- c("lnAUCt","lnAUCinf","lnCmax")
  AUCt.stderr <- attr(lm.lnauct$fixDF,"varFixFact")[2,2]
  AUCt.intervals <- intervals(lm.lnauct,0.90,"fixed")
  AUCt.var <- VarCorr(lm.lnauct)
  AUCinf.stderr <- attr(lm.lnaucinf$fixDF,"varFixFact")[2,2]
  AUCinf.intervals <- intervals(lm.lnaucinf,0.90,"fixed")
  AUCinf.var <- VarCorr(lm.lnaucinf)
  Cmax.stderr <- attr(lm.lncmax$fixDF,"varFixFact")[2,2]
  Cmax.intervals <- intervals(lm.lncmax,0.90,"fixed")
  Cmax.var <- VarCorr(lm.lncmax)
  
  Stderr <- c(AUCt.stderr, AUCinf.stderr, Cmax.stderr)
  LowerCI <- c(AUCt.intervals$fixed[2,1], AUCinf.intervals$fixed[2,1],Cmax.intervals$fixed[2,1])
  Ratio <- c(AUCt.intervals$fixed[2,2], AUCinf.intervals$fixed[2,2],Cmax.intervals$fixed[2,2])
  UpperCI <- c(AUCt.intervals$fixed[2,3], AUCinf.intervals$fixed[2,3],Cmax.intervals$fixed[2,3])
  
  ANOVA_lnResults <- data.frame(Parameter, Stderr, Ratio, LowerCI, UpperCI)
  
  Parameter <- c("AUCt","AUCinf","Cmax")
  IntraCV <- c(100*((exp(as.numeric(AUCt.var[2]))-1)^0.5),100*((exp(as.numeric(AUCinf.var[2]))-1)^0.5),100*((exp(as.numeric(Cmax.var[2]))-1)^0.5))
  InterCV <- c(100*((exp(as.numeric(AUCt.var[1]))-1)^0.5),100*((exp(as.numeric(AUCinf.var[1]))-1)^0.5),100*((exp(as.numeric(Cmax.var[1]))-1)^0.5))
  LowerCI <- 100*exp(LowerCI)
  Ratio <- 100*exp(Ratio)
  UpperCI <- 100*exp(UpperCI)
  success<- (LowerCI > 80) & (UpperCI < 125)
  pVal_TRT <- c(summary(lm.lnauct)$tTable[2,5],summary(lm.lnaucinf)$tTable[2,5],summary(lm.lncmax)$tTable[2,5])
  pVal_PER <- c(summary(lm.lnauct)$tTable[3,5],summary(lm.lnaucinf)$tTable[3,5],summary(lm.lncmax)$tTable[3,5])
  pVal_SEQ <- c(summary(lm.lnauct)$tTable[4,5],summary(lm.lnaucinf)$tTable[4,5],summary(lm.lncmax)$tTable[4,5])
  
  #Calculate Power
  n_auct <- lm.lnauct$dims$N/2
  tau1_auct <- sqrt(n_auct)*log(Ratio[1]/80)/(sqrt(2*as.numeric(AUCt.var[2])))
  tau2_auct <- sqrt(n_auct)*log(Ratio[1]/125)/(sqrt(2*as.numeric(AUCt.var[2])))
  Power_auct <- 100*(pt(tau1_auct-qt(0.95,n_auct-2),df=n_auct-2) - pt(tau2_auct+qt(0.95,n_auct-2),df=n_auct-2))
  
  n_aucinf <- lm.lnaucinf$dims$N/2
  tau1_aucinf <- sqrt(n_aucinf)*log(Ratio[2]/80)/(sqrt(2*as.numeric(AUCinf.var[2])))
  tau2_aucinf <- sqrt(n_aucinf)*log(Ratio[2]/125)/(sqrt(2*as.numeric(AUCinf.var[2])))
  Power_aucinf <- 100*(pt(tau1_aucinf-qt(0.95,n_aucinf-2),df=n_aucinf-2) - pt(tau2_aucinf+qt(0.95,n_aucinf-2),df=n_aucinf-2))
  
  n_cmax <- lm.lncmax$dims$N/2
  tau1_cmax <- sqrt(n_cmax)*log(Ratio[3]/80)/(sqrt(2*as.numeric(Cmax.var[2])))
  tau2_cmax <- sqrt(n_cmax)*log(Ratio[3]/125)/(sqrt(2*as.numeric(Cmax.var[2])))
  Power_cmax <- 100*(pt(tau1_cmax-qt(0.95,n_cmax-2),df=n_cmax-2) - pt(tau2_cmax+qt(0.95,n_cmax-2),df=n_cmax-2))
  
  Power <- c(Power_auct, Power_aucinf, Power_cmax)
  
  ANOVA_Results <- data.frame(Parameter, IntraCV, InterCV, Ratio, LowerCI, UpperCI, pVal_TRT, pVal_PER, pVal_SEQ, Power,success)
  
  return(ANOVA_Results)
}
bootstrap_ABE <- function(data, n_subjects, n_bootstrap, 
                          output_file = "bootstrap_results.csv",
                          output_file_sampling = "bootstrap_sampling.csv") {
  
  results <- data.frame()
  sampling_info <-data.frame()
  success_count <- 0
  
  for (i in 1:n_bootstrap) {
    # Resample subjects with replacement
    sampled_subjects <- sample(unique(data$SUBJ), size = n_subjects, replace = TRUE)
    # Extract records for sampled subjects
    sampled_data<-data.frame()# rest the data at the beginning of each loop
    
    for ( j in 1:n_subjects) {
      
      dataj <- data[data$SUBJ %in% sampled_subjects[j], ]%>%
        mutate(SUBJ = j) 
      sampled_data <- rbind(sampled_data,dataj)
      
    }
    # Run ABE analysis
    # abe_results <- abe_analysis(sampled_data)
    # itration_results <- data.frame(
    #  Iteration = rep(i,3),
    # abe_results)
    
    aberesult <- try(abe_analysis(sampled_data), silent = TRUE)
    if (class(aberesult) == "try-error") {
      # 如果函数运行失败，记录错误信息
      print(paste("Error occurred in iteration", i))
      next
    } 
    
    else 
    {
      # 如果函数运行成功，增加成功次数计数器
      success_count <- success_count + 1
      itration_results <- data.frame(
        Iteration = rep(i,3),
        abe_analysis(sampled_data))
    }
    
    
    
    iteration_subjects <- data.frame(
      Iteration = rep(i,length(sampled_subjects)),
      Subjects = sampled_subjects
    )
    
    
    # Save results,在循环外每次更新结果数据
    sampling_info <- rbind(sampling_info, iteration_subjects)
    results <- rbind(results,itration_results)
  }
  #write.csv(results,file= output_file, row.names = FALSE)
  #write.csv(sampling_info,file= output_file_sampling, row.names = FALSE)
  
  print(paste("Bootstrap success count is:", success_count))
  return(results)
  
}

shinyApp(ui = ui, server = server)
