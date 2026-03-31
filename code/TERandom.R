#### te.pipeline 2025-10-19
#### TE Random
install.packages("ranger")
install.packages("mlr3")
remove.packages("nanonext") 
require(devtools)
install.packages("nanonext", version = "1.7.2")
install.packages("bbotk", version = "1.8.0")
install.packages("rlang", version = "1.1.4")
install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, mlr3verse, data.table)
pacman::p_load(mlr3verse)
install.packages("vip")
##
library(ranger)
library(tidyverse) 
library(mlr3)
library(mlr3verse)
if(T){
  group <- factor(c(rep("B3",20),rep("B1",20),rep("B2",20),rep("B4",20),
                    rep("H3",20),rep("H1",20),rep("H2",20),rep("H4",20),
                    rep("L3",20),rep("L1",20),rep("L2",20),rep("L4",20),
                    rep("M3",20),rep("M1",20),rep("M2",20),rep("M4",20),
                    rep("S3",20),rep("S1",20),rep("S2",20),rep("S4",20)))
  group2 <- factor(c(rep('Middle_Aged',20),rep('Adolescent',20),rep('Adult',20),rep('Old',20)))
}
#family.class
tissue ='Brain'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
family.class <- result.k[,c(1,2)]
##brain
tissue ='Brain'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
m <- as.character(group2)
m <- m[-c(61,62)]
if(T){
result.k <- result.k[,-1]
result2 <- result.k %>% group_by(family_id) %>% 
  dplyr::summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k$tissue <- tissue
result.k$group <- m
}
result.b <- result.k #brain
##heart
tissue ='Heart'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
m <- as.character(group2)
m <- m[-c(62)] 
if(T){
result.k <- result.k[,-1]
result2 <- result.k %>% group_by(family_id) %>% 
  dplyr::summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k$tissue <- tissue
result.k$group <- m
}
result.h <- result.k #heart
##lung
tissue ='Lung'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
m <- as.character(group2)
#m <- m[-c(62)] 
if(T){
result.k <- result.k[,-1]
result2 <- result.k %>% group_by(family_id) %>% 
  dplyr::summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k$tissue <- tissue
result.k$group <- m
}
result.l <- result.k #lung
##muscle
tissue ='Muscle'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
m <- as.character(group2)
m <- m[-c(1,4)]
if(T){
result.k <- result.k[,-1]
result2 <- result.k %>% group_by(family_id) %>% 
  dplyr::summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k$tissue <- tissue
result.k$group <- m
}
result.m <- result.k #muscle
##skin
tissue ='Skin'
if(T){
  path <- 'limma_result/'
  rds11 <- paste0(path,tissue,'_family.rds')
}
result.k <- readRDS(rds11)
m <- as.character(group2)
m <- m[-c(23,25,33,35,48,50)]
if(T){
result.k <- result.k[,-1]
result2 <- result.k %>% group_by(family_id) %>% 
  dplyr::summarise(across(everything(), sum,na.rm = TRUE))
result.k <- as.data.frame(result2)
rownames(result.k) <- result.k[,1]
result.k <- result.k[,-1]
result.k <- as.data.frame(t(result.k))
result.k$tissue <- tissue
result.k$group <- m
}
result.s <- result.k #skin
##all
result <- rbind(result.b,result.h)
result <- rbind(result,result.l)
result <- rbind(result,result.m)
result <- rbind(result,result.s)

df = as.data.table(result)
head(df)
saveRDS(df,'random.rds')
df <- readRDS('random.rds')
df <-df[,c('Alu','B4','B2','MIR',
                       'ID','5S-Deu-L2',
                       'L1','L2','CR1','RTE-Bo.B',
                       'Penelope','I-Jockey',
                       'ERVK','ERVL-MaLR',
                       'ERVL','ERV1',
                       'LTR','Gypsy','tissue','group')]
df2 <- df
table(df2$tissue)#Brain  Heart   Lung Muscle   Skin 

df <- df2[which(df2$tissue=='Skin'),]
df <-df[,c('Alu','B4','B2','MIR',
           'ID','5S-Deu-L2',
           'L1','L2','CR1','RTE-Bo.B',
           'Penelope','I-Jockey',
           'ERVK','ERVL-MaLR',
           'ERVL','ERV1',
           'LTR','Gypsy','group')]
# 1. 创建分类任务 (Task)
# 指定目标变量为 "Species"
task = as_task_classif(df, target = "group", id = "iris_task")
# 可视化任务，查看特征分布
autoplot(task)

# 2. 划分训练集和测试集
set.seed(42) # 设置随机种子以保证结果可重现
split = partition(task, ratio = 0.8) # 80% 作为训练集，20% 作为测试集
cat("训练集样本数:", length(split$train), "\n")
cat("测试集样本数:", length(split$test), "\n")

# 3.构建图形化管道（预处理 + 模型）
# **预处理管道 (prep)**
prep = po("fixfactors") %>>%         # 修正因子型特征，确保一致性
  po("removeconstants") %>>% # 移除常量特征
  po("scale")                # 对数值型特征进行缩放/标准化
# **学习器 (rf)：随机森林分类器**
# 随机森林是一个强大的集成学习算法
rf = lrn("classif.ranger", #
         predict_type = "prob",
         importance = "impurity", # 启用特征重要性计算
         # 标记需要调优的超参数，使用 to_tune() 函数定义搜索范围
         num.trees = to_tune(200, 800), 
         mtry = to_tune(p_int(2, 4))) # mtry: 每次分裂随机采样的特征数量
if(F){
learners <- list( 
  "Logistic" = lrn("classif.log_reg", predict_type = "prob"), 
  "决策树" = lrn("classif.rpart", predict_type = "prob"), 
  "随机森林" = lrn("classif.ranger", predict_type = "prob"), 
  "SVM" = lrn("classif.svm", type = "C-classification", predict_type = "prob"), 
  "神经网络" = lrn("classif.nnet", predict_type = "prob"), 
  "XGBoost" = lrn("classif.xgboost", predict_type = "prob"), 
  "朴素贝叶斯" = lrn("classif.naive_bayes", predict_type = "prob"), 
  "KNN" = lrn("classif.kknn", predict_type = "prob") )
rf <- learners[[1]]
}
# **组装图形化学习器**
graph = prep %>>% rf
glrn = as_learner(graph) # 将管道转换为一个可训练的 Learner 对象
# 可视化完整的管道结构
graph$plot(horizontal = TRUE)

# 4.超参数调优与模型训练
# **超参数调优 (at) 规则设定**
at = auto_tuner(
  tuner = tnr("random_search"),          # 调优方法：随机搜索
  learner = glrn,                        # 待调优的图形化学习器
  resampling = rsmp("cv", folds = 5),    # 评估方法：5折交叉验证、3
  measure = msr("classif.mauc_au1u"),    # 优化指标：多类别 AUC
  term_evals = 15)                       # 搜索预算：最多评估 15 次
# **执行训练与调优**
#future::plan("multisession") # 启用并行计算加速
#set.seed(1)
at$train(task, row_ids = split$train) # 在训练集上进行调优和最终模型训练
# 查看调优结果
at$tuning_result

# 5.性能评估与可视化
#A. 特征重要性分析
# 查看最佳超参数和特征重要性
print(at$tuning_result[which.min(at$tuning_result$classif.mauc_au1u),])
# 可视化特征重要性
# 识别对预测鸢尾花物种贡献最大的特征

library(vip)
at$base_learner() |> 
  vip::vip(aesthetics = list(fill = "#CD6868"))

#B. 测试集预测与混淆矩阵
# 在测试集上进行预测
pred = at$predict(task, row_ids = split$test)
# 混淆矩阵 (Confusion Matrix)
# 混淆矩阵清晰展示了模型对每个物种的分类准确度，以及错误是如何发生的。
# 对角线上的数值是正确分类的样本数。
cm = pred$confusion
cm

# C. 评估指标计算
# 评估指标
msrs(c("classif.acc",         # 准确率 (Accuracy)
       "classif.mauc_au1u",   # 多类别 AUC (越接近 1 越好)
       "classif.mbrier")) |>  # 多类别 Brier 分数 (越接近 0 越好)
  pred$score()

# D. ROC 曲线可视化
# ROC 曲线（多分类）
library(yardstick) # 用于多分类评估和可视化
library(ggplot2)
library(data.table)
if(F){
as.data.table(pred) |>
  roc_curve(truth = truth, 
            prob.Adolescent, prob.Adult, prob.Middle_Aged, prob.Old) |> # 传入真实值和每个类别的概率
  autoplot() + 
  labs(title = "ROC") # 曲线越靠近左上角，性能越好
autoplot(pred,type='roc')
}
library(data.table)
library(pROC)
library(ggplot2)

# 假设 pred 是包含预测概率和真实标签的数据框
# 转换为 data.table 格式
df <- as.data.table(pred)

# 绘制 ROC 曲线
# 方法3: 计算每个类别AUC（使用二分类方式）
df_0vsrest <- df %>%
  mutate(truth_binary = ifelse(truth == 'Adolescent',0,1))
df_1vsrest <- df %>%
  mutate(truth_binary = ifelse(truth == 'Adult',0,1))
df_2vsrest <- df %>%
  mutate(truth_binary = ifelse(truth == 'Middle_Aged',0,1))
df_3vsrest <- df %>%
  mutate(truth_binary = ifelse(truth == 'Old',0,1))
# 计算每个类别的AUC
df_0vsrest$truth_binary <- as.factor(df_0vsrest$truth_binary)
df_1vsrest$truth_binary <- as.factor(df_1vsrest$truth_binary)
df_2vsrest$truth_binary <- as.factor(df_2vsrest$truth_binary)
df_3vsrest$truth_binary <- as.factor(df_3vsrest$truth_binary)
auc_class0 <- df_0vsrest %>%
  roc_auc(truth_binary, prob.Adolescent, estimator = "binary")
auc_class1 <- df_1vsrest %>%
  roc_auc(truth_binary, prob.Adult, estimator = "binary")
auc_class2 <- df_2vsrest %>%
  roc_auc(truth_binary, prob.Middle_Aged,estimator = "binary")
auc_class3 <- df_3vsrest %>%
  roc_auc(truth_binary, prob.Old, estimator = "binary")
auclist <- c(paste0('AUC (','Adolescent','): ',round(auc_class0$.estimate,5)),
             paste0('AUC (','Adult','): ',round(auc_class1$.estimate,5)),
             paste0('AUC (','Middle_Aged','): ',round(auc_class2$.estimate,5)),
             paste0('AUC (','Old','): ',round(auc_class3$.estimate,5))
             )

# 计算多分类 ROC 曲线
roc_data <- df %>% 
  roc_curve(truth = truth, 
            prob.Adolescent, prob.Adult, prob.Middle_Aged, prob.Old)
colnames(roc_data) <- c("group","threshold","specificity","sensitivity")

# 绘制ROC曲线图
theme_dlqpect <- function(...){
  theme(
    panel.grid = element_blank(),
    axis.text =  element_text(size = 9, family = 'Times New Roman'),
    axis.title = element_text(size = 10, family = 'Times New Roman'),
    plot.title = element_text(hjust = 0.5, vjust = 1,size = 10,
                              margin = margin(t = 10),face = "bold",
                              family = 'Times New Roman'),
    legend.text = element_text(size = 9,,family = "Times New Roman"),
    legend.title = element_blank() )
}
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = group)) +
  geom_path(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  labs(x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       title = "ROC Curves") +
  annotate(geom = "text", x= 0.75, y=  0.1, color='#cd3232',
           label=auclist[1],   size=3,family = "Times New Roman")+
  annotate(geom = "text", x= 0.75, y=  0.2, color='#e59f00',
           label=auclist[2],   size=3,family = "Times New Roman")+
  annotate(geom = "text", x= 0.75, y=  0.3, color='#cd7ba9',
           label=auclist[3],   size=3,family = "Times New Roman")+
  annotate(geom = "text", x= 0.75, y=  0.4, color='#019e74',
           label=auclist[4],   size=3,family = "Times New Roman")+
  theme_minimal() +theme_bw()+theme_dlqpect()+
  scale_color_manual(" ",values =c('#cd3232','#e59f00','#cd7ba9','#019e74'))+
  theme(legend.position = "bottom")


