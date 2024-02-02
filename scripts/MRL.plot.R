#MRL.plot.R

#Examines transcript levels of genes-of-interest (e.g MRL) based on publicly available mRNA-seq datasets

dirRoot = '/Users/michael.nodine/Documents/manuscripts/primary/Honglei_MRL/github/'
dataRoot = paste(dirRoot,'resources/',sep="")
printRoot = paste(dirRoot,'graphs/',sep="")

get.expression.global <- function() {
  df=read.delim(paste(dataRoot,"mRNA_TPM.tsv",sep=""), sep='\t',
                row.names=1, header=TRUE, stringsAsFactors=FALSE,strip.white=TRUE)
  return(df)
}

plot.targets.ext.pdf.global <- function(sRNA,tar.v,y.max) {
  print(y.max)
  print(sRNA)
  print(tar.v)
  
  get.exp.sub <- function(all.df) {
    df.tar=subset(all.df, rownames(all.df) %in% tar.v ==TRUE)
    return(df.tar)
  }
  
  #General development
  df.tar.global=get.exp.sub(get.expression.global())
  
  for (tar in tar.v) {
    print(tar)
    
    
    global.list=list("2cell"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[1:2]),
                     "8cell"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[3:4]),
                     "pg"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[5:7]),
                     "32cell"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[8:9]),
                     "gl"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[10:12]),
                     "eh"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[13:15]),
                     "lh"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[16:18]),
                     "et"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[19:21]),
                     "lt"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[22:24]),
                     "emb_7dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[25:27]),
                     "emb_8dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[28:30]),
                     "bc"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[31:33]),
                     "emb_10dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[34:36]),
                     "mg"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[37:39]),
                     "emb_12dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[40:42]),
                     "emb_13dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[43:45]),
                     "emb_15dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[46:48]),
                     "emb_17dap"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[49:50]),
                     "SE_5D"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[51]),
                     "SE_10D"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[52]),
                     "SE_15D"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[53]),
                     "callus_neg"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[54:55]),
                     "callus_lec1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[56:57]),
                     "seed_dai0"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[58:59]),
                     "seed_dai1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[60:61]),
                     "seed_dai2"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[62:63]),
                     "seed_dai3"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[64:65]),
                     "cot_1dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[66:67]),
                     "hycot_1dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[68:69]),
                     "root_1dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[70:71]),
                     "root_7dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[72:73]),
                     "root.apex_7dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[74:75]),
                     "sam_1dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[76:77]),
                     "sam_7dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[78]),
                     "sam_8dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[79]),
                     "sam_9dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[80:81]),
                     "sam_10dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[82:83]),
                     "sam_11dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[84:85]),
                     "sam_12dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[86:87]),
                     "sam_13dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[88:89]),
                     "sam_14dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[90]),
                     "sam_15dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[91]),
                     "sam_16dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[92]),
                     "petiole_7dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[93:94]),
                     "petiole_9dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[95:96]),
                     "petiole_12dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[97:98]),
                     "petiole_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[99:100]),
                     "petiole_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[101:102]),
                     "blade_7dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[103:104]),
                     "blade_9dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[105:106]),
                     "blade_12dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[107:108]),
                     "blade_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[109:110]),
                     "leafvein_12dag"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[111:112]),
                     "leafvein_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[113:114]),
                     "leafvein_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[115:116]),
                     "leaf_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[117:118]),
                     "leaf_SL"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[119:120]),
                     "internode_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[121:122]),
                     "internode_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[123:124]),
                     "pedicel_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[125:126]),
                     "floraxis_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[127:128]),
                     "inflor_mat"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[129:130]),
                     "fb_SL"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[131:132]),
                     "fb_MK"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[133:135]),
                     "flower_st8"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[136:137]),
                     "flower_st9"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[138:139]),
                     "flower_st10"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[140:141]),
                     "flower_st11"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[142:143]),
                     "flower_st12.1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[144:145]),
                     "flower_st12.2"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[146:147]),
                     "flower_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[148:149]),
                     "flower_st14"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[150:151]),
                     "flower_st15"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[152:153]),
                     "sepal_st9"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[154:155]),
                     "sepal_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[156:157]),
                     "petal_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[158:159]),
                     "filament_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[160:161]),
                     "anther_st9"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[162:163]),
                     "anther_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[164:165]),
                     "anther_st15"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[166:167]),
                     "stigma_st11"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[168:169]),
                     "pistil_st9"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[170:171]),
                     "pistil_st13"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[172:173]),
                     "silique1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[174:175]),
                     "silique2"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[176:177]),
                     "silique3"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[178:179]),
                     "silique4"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[180:181]),
                     "silique_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[182:183]),
                     "carpel"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[184:185]),
                     "silpod1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[186:187]),
                     "silpod2"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[188:189]),
                     "silpod3"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[190:191]),
                     "silpod4"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[192:193]),
                     "silpod_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[194:195]),
                     "ovule"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[196:197]),
                     "seed1"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[198:199]),
                     "seed2"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[200:201]),
                     "seed3"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[202:203]),
                     "seed4"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[204:205]),
                     "seed5"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[206:207]),
                     "seed6"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[208:209]),
                     "seed7"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[210:211]),
                     "seed8"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[212:213]),
                     "seed9"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[214:215]),
                     "seed_sen"=subset(df.tar.global,rownames(df.tar.global) %in% tar == TRUE, select=names(df.tar.global)[216:217]))
    
    add.arrows <- function(exp.list,diff) {
      
      for (i in c(1:length(names(exp.list)))) {
        arrows(i-diff, mean(as.numeric(exp.list[[i]])), i+diff, mean(as.numeric(exp.list[[i]])), code=0, angle=90, lwd=2)
      }
    }
    
    plotFile=paste(printRoot,sRNA,".",tar,".",y.max,".all.pdf",sep="")
    print(plotFile)
    pdf(plotFile, useDingbats=FALSE, width=25, height=5)
    par(mar=c(7,4,4,2) + 0.1)
    #global.list
    if (length(unlist(global.list)) > 0) {
      if (y.max == "auto") {
        print("auto ymax")
        ymax=max(unlist(global.list))
        print(ymax)
        stripchart(global.list,ylab="mRNA levels (TPM)", vertical=TRUE, method='jitter', pch=19, las=2, col="black", 
                   xlab="", main=paste(tar," global",sep=""), ylim=c(0,ymax))
        add.arrows(global.list,0.2)
        abline(h=1,lty=2)
      }
      if (y.max != "auto") {
        ymax=y.max
        print("set ymax")
        print(ymax)
        stripchart(global.list,ylab="mRNA levels (TPM)", vertical=TRUE, method='jitter', pch=19, las=2, col="black", 
                   xlab="", main=paste(tar," global",sep=""), ylim=c(0,ymax))
        add.arrows(global.list,0.2)
        abline(h=1,lty=2)
      }
    }
    dev.off()
  }
}

###########################################################################################################################################
plot.targets.ext.pdf.global("MRL",c("AT5G14050"),100)
plot.targets.ext.pdf.global("AT1G01830",c("AT1G01830"),100)
plot.targets.ext.pdf.global("MRL",c("AT5G14050"),"auto")
plot.targets.ext.pdf.global("AT1G01830",c("AT1G01830"),"auto")