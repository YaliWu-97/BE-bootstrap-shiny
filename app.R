library(shiny)
library(nlme)
library(dplyr)
library(tidyr)
library(DT)
select <- dplyr::select
filter <- dplyr::filter
#functions used for x2way crossover####
abe_analysis <- function(data){
  ABEexample1 <- data
  
  
  p <- 3 # rounding for Cmax and AUCs
  pci <- 2 # rounding for ratio and 90%CI (use 1 for TPD, 2 for EMA, FDA)
  
  
  factors       <- c("SUBJ", "GRP", "PRD", "TRT")
  ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for lme()
  
  lm.lnauct <- lme(log(AUClast)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  lm.lnaucinf <- lme(log(AUCinf)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  lm.lncmax <- lme(log(Cmax)~TRT+PRD+GRP, data=ABEexample1, random=~1|SUBJ,na.action=na.exclude)
  
  
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
  return(list(results = results, sampling_info = sampling_info))
  
}
#functions used for x3x4 crossover####
#data, data2是未经censor的数据
data.censored <- function(data,para,n){
  # 找出有缺失值的个体
  subj_with_na <- unique(data$SUBJ[!complete.cases(data[, para])])
  
  # 删除这些个体的所有数据
  data_clean <- data[!data$SUBJ %in% subj_with_na, ]
  
  #检查每个受试者的记录数
  record_count <- table(data$SUBJ)
  
  #找出记录数不足的受试者
  incomplete_subj <- names(record_count[record_count < n])  # 假设完整记录数为4
  
  #删除这些受试者的所有记录
  data_clean <- data[!data$SUBJ %in% incomplete_subj, ]
  return(data_clean)
}
scavbe_f <- function(data,n_crossover,para){
  crossover_test<- data%>%
    select(SUBJ,GRP,PRD,TRT,!!sym(para))%>%#dplyrselect的xin
    filter(TRT=="T")%>%
    mutate(LOGpara=log(!!sym(para)))%>%##动态引用参数
    select(SUBJ,GRP,PRD,TRT,LOGpara)
  
  crossover_ref<- data%>%
    select(SUBJ,GRP,PRD,TRT,!!sym(para))%>%
    filter(TRT=="R")%>%
    mutate(LOGpara=log(!!sym(para)))%>%#    
    select(SUBJ,GRP,PRD,TRT,LOGpara) 
  if (n_crossover == 3){
    ref_mean <- crossover_ref%>%
      group_by(SUBJ)%>%
      summarise(ref_m = mean(LOGpara))
    dlat <- crossover_ref%>%
      group_by(SUBJ)%>%
      mutate (dlat = diff(LOGpara))%>%
      select(SUBJ,dlat)%>%
      distinct(SUBJ,dlat)
    
    scavbe <- merge(ref_mean,dlat, by=c("SUBJ"))%>%
      merge(crossover_test,by=c("SUBJ"))%>%
      mutate(ilat = as.numeric(LOGpara) - as.numeric(ref_m))
    
    # 待有三周期数据集时进行补充
  }else if (n_crossover ==4){
    lat1t<-subset(crossover_test, PRD == 1 | PRD == 2)
    lat2t<-subset(crossover_test, PRD == 3 | PRD == 4)
    lat1r<-subset(crossover_ref, PRD == 1 | PRD == 2)
    lat2r<-subset(crossover_ref, PRD == 3 | PRD == 4)
    
    colnames(lat1t)[colnames(lat1t) == "LOGpara"] <- "LOGparaT1"
    colnames(lat2t)[colnames(lat2t) == "LOGpara"] <- "LOGparaT2"
    colnames(lat1r)[colnames(lat1r) == "LOGpara"] <- "LOGparaR1"
    colnames(lat2r)[colnames(lat2r) == "LOGpara"] <- "LOGparaR2"
    
    scavbe<-merge(lat1t,lat2t,by=c("SUBJ","GRP"))%>%
      merge(lat1r,by=c("SUBJ","GRP"))%>%
      merge(lat2r,by=c("SUBJ","GRP"),suffixes = c(".left", ".right"))%>%
      #防止列明重复总是有warning，循环时比较麻烦
      select(SUBJ,GRP,LOGparaT1,LOGparaT2,LOGparaR1,LOGparaR2)%>%
      mutate(ilat=0.5*(LOGparaT1+LOGparaT2-LOGparaR1-LOGparaR2))%>%
      mutate(dlat=LOGparaR1-LOGparaR2)
    factors       <- c("SUBJ", "GRP")
    scavbe[factors] <- lapply(scavbe[factors], factor) # factorize for lme()
    
  } else{
    print("You can only enter 3 or 4 as crossover number in this function")
  }
  #    return(crossover_test)
  #    return(crossover_ref)
  return(scavbe)
}
abe_highV<-function(data,para){
  formula0 <- as.formula(paste0(para, "~ TRT | SUBJ"))
  data<- groupedData(formula0, data = data)
  data<-data%>%
    select(SUBJ,GRP,PRD,TRT,!!sym(para))%>%
    mutate(para = !!sym(para))%>%
    select(SUBJ,GRP,PRD,TRT,para)
  factors       <- c("SUBJ", "GRP", "PRD", "TRT")
  data[factors] <- lapply(data[factors], factor) # factorize for lme()
  
  lm.lnpara.ABE <- lme(log(para) ~ TRT + GRP + PRD, 
                       random = list(SUBJ = pdBlocked(list(pdIdent(~1),
                                                           pdIdent(~TRT-1)))) ,
                       data = data,
                       #weights = varIdent(form = ~ 1 | TRT),
                       correlation = corSymm( value = 0.3, form = ~ 1|SUBJ/TRT),#若出现不收敛的情况可以在这里更改初始值
                       method = "REML"
  )
  
  intervals(lm.lnpara.ABE,0.90,"fixed")
  Cmax.intervals <- intervals(lm.lnpara.ABE,0.90,"fixed")
  LowerCI <- c(Cmax.intervals$fixed[2,1])
  UpperCI <- c(Cmax.intervals$fixed[2,3])
  LowerCI <- 100* exp(LowerCI)
  UpperCI <- 100* exp(UpperCI)
  CI_logform<-c(LowerCI,UpperCI)
  return(CI_logform)
}
RSABE<- function(data2,n_crossover,para){
  
  data<- data.censored(data2,para,n_crossover)
  scavbe<-scavbe_f(data,n_crossover,para)
  lm_model <- lm(dlat ~ GRP, data = scavbe)
  ref_std2<-summary(lm_model)$sigma^2#残差标准差，它是残差平方和的平方根。
  s2wr<-ref_std2/2# 根据指导原则，平方后除以2得到s2wr
  dfd <- as.numeric(lm_model$df.residual)
  #calcultate the pointest,x,boundx from scavbe dataset
  #using ilat 
  
  lm.lnpara<-lm(ilat~GRP, 
                data=scavbe, 
                na.action=na.exclude)
  #para.stderr <- attr(lm.lnpara$fixDF,"varFixFact")[2,2]#固定效应序列的std
  summary_model <- summary(lm.lnpara)
  se.fit<-predict(lm.lnpara,se.fit = TRUE)$se.fit
  len<-length(se.fit)
  var<-sum(se.fit^2)/len
  if (n_crossover == 3){#if crossover=3,var/3,if crossover=4,var/2
    para.stderr <- sqrt(var/3)
  }
  else if (n_crossover==4){
    para.stderr <- sqrt(var/2) 
  }
  
  
  
  #判断是三周期还是四周期，此处的点估计值不同
  if(n_crossover == 3){
    estimate <- lm.lnpara$coefficients[1]+
      lm.lnpara$coefficients[2]*0.33333 +
      lm.lnpara$coefficients[3]*0.33333
    
  }else if(n_crossover == 4){
    
    estimate <- lm.lnpara$coefficients[1]+lm.lnpara$coefficients[2]*0.5
    #四周期为两个序列均值，若为3周期则为0.3333*对应序列
    
  }
  
  pointest <-  exp(estimate)
  x=estimate**2- para.stderr**2;
  LowerCI <- estimate - qt(0.95, dfd) * para.stderr
  UpperCI <- estimate + qt(0.95, dfd) * para.stderr
  
  boundx=(max((abs(LowerCI)),(abs(UpperCI))))**2;
  
  #置信区间的计算及判定
  theta=((log(1.25))/0.25)**2;
  y=-theta*s2wr;
  boundy=y*dfd/qchisq(0.95,dfd)#
  SWR=sqrt(s2wr);
  criticalbound=(x+y)+sqrt(((boundx-x)**2)+((boundy-y)**2))
  
  LowerCI<-100*exp(LowerCI)
  UpperCI<-100*exp(UpperCI)
  
  
  if(SWR<0.294){
    LowerCI<-abe_highV(data2,para)[1]
    UpperCI<-abe_highV(data2,para)[2]
  }
  #在此处判断给出的CI是否需要由ABE给出
  
  result<-data.frame(para, LowerCI,UpperCI,SWR,pointest,criticalbound)
  
  return(result)
}

bootstrap_RSABE <- function(data, n_subjects, n_bootstrap, n_crossover, params_list,
                            output_file = "bootstrap_results.csv",
                            output_file_sampling = "bootstrap_sampling.csv") {
  
  # 初始化结果列表
  all_results <- data.frame()
  sampling_info <- data.frame()
  
  # 初始化每个参数的成功次数
  success_counts <- setNames(rep(0, length(params_list)), params_list)
  
  # 遍历 Bootstrap 次数
  for (i in 1:n_bootstrap) {
    # Resample subjects with replacement
    sampled_subjects <- sample(unique(data$SUBJ), size = n_subjects, replace = TRUE)
    
    # Extract records for sampled subjects
    sampled_data <- data.frame()  # 重置数据
    
    for (j in 1:n_subjects) {
      dataj <- data[data$SUBJ %in% sampled_subjects[j], ] %>%
        mutate(SUBJ = j)
      sampled_data <- rbind(sampled_data, dataj)
    }
    
    # 初始化当前 Bootstrap 迭代的结果
    iteration_results <- data.frame()
    
    # 遍历参数列表
    for (param in params_list) {
      # 对每个参数调用 RSABE 函数
      rsaberesult <- tryCatch({
        RSABE(sampled_data, n_crossover, param)
      }, error = function(e) {
        print(paste("Error occurred in iteration", i, "for parameter", param))
        return(NULL)  # 如果出现错误，返回 NULL
      })
      
      if (!is.null(rsaberesult)) {
        # 添加参数名称列和迭代次数列
        rsaberesult$Parameter <- param
        rsaberesult$Iteration <- i
        
        # 判断 BE 实验是否成功
        rsaberesult <- rsaberesult %>%
          mutate(
            success = case_when(
              SWR < 0.294 & LowerCI >= 80 & UpperCI <= 125 ~ TRUE,
              SWR >= 0.294 & pointest >= 0.8 & pointest <= 1.25 & criticalbound < 0 ~ TRUE,
              TRUE ~ FALSE
            )
          )
        
        # 更新当前参数的成功次数
        if (rsaberesult$success) {
          success_counts[param] <- success_counts[param] + 1
        }
        
        # 添加当前 Bootstrap 迭代的结果
        iteration_results <- rbind(iteration_results, rsaberesult)
      }
    }
    
    # 合并当前 Bootstrap 迭代的结果
    if (nrow(iteration_results) > 0) {
      all_results <- rbind(all_results, iteration_results)
      
      # 保存抽样信息
      iteration_subjects <- data.frame(
        Iteration = rep(i, length(sampled_subjects)),
        Subjects = sampled_subjects
      )
      sampling_info <- rbind(sampling_info, iteration_subjects)
    }
  }
  
  # 打印每个参数的成功次数
  for (param in params_list) {
    print(paste("Bootstrap success count for", param, "is:", success_counts[param]))
  }
  
  # 返回结果
  return(list(results = all_results, sampling_info = sampling_info))
}










#ui setup ####

ui <- fluidPage(
  titlePanel("BE analysis tools"),  # 标题
  sidebarLayout(
    sidebarPanel(
      # 选择试验类型
      radioButtons("trial_type", "choose your trail type",
                   choices = c("2×2 crossover study", "Partial or fully replicated crossover"),
                   selected = "2×2 crossover study"),  # 默认选择
      # 示例表格
      h5("Example Data Format"),
      img(src = "example.png", width = "100%"),

      # 动态 UI：根据试验类型显示不同的输入控件
      uiOutput("dynamic_ui")
    ),
    mainPanel(
      # 动态 UI：根据试验类型显示不同的结果
      uiOutput("dynamic_output")
    )
  )
)




#server setup####
# 动态生成输入控件
server <- function(input, output) {

  output$dynamic_ui <- renderUI({
    if (input$trial_type == "2×2 crossover study") {
      # 2×2 交叉试验的输入控件
      tagList(
        fileInput("file", "upload CSV file", accept = ".csv"),  # 文件上传

        numericInput("n_subjects", "number of subjects in each sampling (n_subjects):", value = 10),  # 输入受试者数量
        numericInput("n_bootstrap", "number of Bootstrap(n_bootstrap):", value = 100),  # 输入 Bootstrap 次数
        actionButton("run_abe", "ABE analysis"),  # ABE 分析按钮
        actionButton("run_bootstrap", "run Bootstrap"),  # Bootstrap 分析按钮
        downloadButton("downloadResults", "Download Results"),
        downloadButton("downloadSamplingInfo", "Download Sampling Info")
      )
    } else {
      # 完全重复/部分重复试验的输入控件
      tagList(
        fileInput("file", "upload CSV file", accept = ".csv"),  # 文件上传
        numericInput("n_subjects", "number of subjects in each sampling (n_subjects):", value = 10),  # 输入受试者数量
        numericInput("n_bootstrap", "number of Bootstrap (n_bootstrap):", value = 100),  # 输入 Bootstrap 次数
        numericInput("n_crossover", "enter the period(3 or 4)", value = 3, min = 3, max = 4, step = 1),  # 输入周期数
        textInput("para", "PK parameters(separated by commas).", value = "Cmax, AUClast, AUCINF"),
        actionButton("run_RSABE", "run ABE/RSABE"),  # RSABE 分析按钮
        actionButton("run_bootstrap_RSABE", "run Bootstrap"),  # Bootstrap 分析按钮
        downloadButton("downloadResults", "Download Results"),
        downloadButton("downloadSamplingInfo", "Download Sampling Info")
      )
    }
  })
  
  
  # 动态生成输出
  output$dynamic_output <- renderUI({
    if (input$trial_type == "2×2 crossover study") {
      # 2×2 交叉试验的输出
      tagList(
        h3("ABE results"),  # ABE 结果标题
        tableOutput("abe_result"),  # ABE 结果显示
        h3("Bootstrap results"),  # Bootstrap 结果标题
        tableOutput("bootstrap_result")  # Bootstrap 结果显示
      )
    } else {
      # 完全重复/部分重复试验的输出
      tagList(
        h3("Results of partially/fully replicated study"),  # 结果标题
        tableOutput("RSABE_results"), # 结果显示
        h3("Bootstrap results"),  # 结果标题
        tableOutput("bootstrap_RSABE")  # 结果显示
      )
    }
  })

  bootstrap_original <- reactiveValues(results = NULL, sampling_info = NULL)
  #监听 x2way abe bottom
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
  
  # 监听 x2way Bootstrap 分析按钮
  observeEvent(input$run_bootstrap, {
    req(input$file)  # 确保文件已上传
    
    # 读取上传的 CSV 文件
    data <- read.csv(input$file$datapath)
    
    # 调用 Bootstrap 分析函数
    test_result <- bootstrap_ABE(data, input$n_subjects, input$n_bootstrap)
    
    bootstrap_original$results <- test_result$results
    bootstrap_original$sampling_info <- test_result$sampling_info
    
    bootstrap_result<- test_result$results %>%
      group_by(Parameter) %>%
      summarise(Successcount = sum(success),
                Totalcount = n(),
                success_rate=(Successcount/Totalcount))
    # 显示 Bootstrap 分析结果
    output$bootstrap_result <- renderTable({
      bootstrap_result
    })
    
     # 下载 results CSV 文件
    output$downloadResults <- downloadHandler(
      filename = function() {
        paste("results-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(bootstrap_original$results, file, row.names = FALSE)
      }
    )
    
    # 下载 sampling_info CSV 文件
    output$downloadSamplingInfo <- downloadHandler(
      filename = function() {
        paste("sampling_info-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(bootstrap_original$sampling_info, file, row.names = FALSE)
      }
    )
    
    
    
    
  })
  
  #监听 x3x4 way RSABE bottom
  
  observeEvent(input$run_RSABE, {
    req(input$file)  # 确保文件已上传
    
    # 读取上传的 CSV 文件
    data <- read.csv(input$file$datapath)
    
    # 解析用户输入的药动学参数
    params <- input$para
    params_list <- strsplit(params, ",")[[1]]  # 按逗号分隔
    params_list <- trimws(params_list)  # 去除每个参数前后的空格
    
    # 初始化结果列表
    all_results <- data.frame()
    
    # 遍历参数列表，分别调用 RSABE 函数
    for (param in params_list) {
      result <- RSABE(data, input$n_crossover, param)
      print(result)
      all_results <- bind_rows(all_results, result)
    }
    
  
    # 显示 RSABE 分析结果
    output$RSABE_results <- renderTable({
      all_results %>%
        mutate(across(where(is.numeric), ~ format(., nsmall = 4)))
    })
    
  })
  
  
  
  # 监听 x3x4way Bootstrap 分析按钮
  observeEvent(input$run_bootstrap_RSABE, {
    req(input$file)  # 确保文件已上传
    
    # 读取上传的 CSV 文件
    data <- read.csv(input$file$datapath)
     # 解析用户输入的药动学参数
    params <- input$para
    params_list <- strsplit(params, ",")[[1]]  # 按逗号分隔
    params_list <- trimws(params_list)  # 去除每个参数前后的空格
    
    # 初始化结果列表
    all_results <- data.frame()
    
    # 调用 Bootstrap 分析函数
    test_result <- bootstrap_RSABE(data, input$n_subjects, input$n_bootstrap,
                                   input$n_crossover,params_list)
    
    bootstrap_original$results <- test_result$results
    bootstrap_original$sampling_info <- test_result$sampling_info
    
    bootstrap_result<- test_result$results %>%
      group_by(para) %>%
      summarise(Successcount = sum(success),
                Totalcount = n(),
                success_rate=(Successcount/Totalcount))
    # 显示 Bootstrap 分析结果
    output$bootstrap_RSABE <- renderTable({
      bootstrap_result%>%
        mutate(across(where(is.numeric), ~ format(., nsmall = 4)))
    })
    
    # 下载 results CSV 文件
    output$downloadResults <- downloadHandler(
      filename = function() {
        paste("results-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(bootstrap_original$results, file, row.names = FALSE)
      }
    )
    
    # 下载 sampling_info CSV 文件
    output$downloadSamplingInfo <- downloadHandler(
      filename = function() {
        paste("sampling_info-", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(bootstrap_original$sampling_info, file, row.names = FALSE)
      }
    )
    
    
    
    
  })
  
  
}







#shinyApp####
shinyApp(ui = ui, server = server)
