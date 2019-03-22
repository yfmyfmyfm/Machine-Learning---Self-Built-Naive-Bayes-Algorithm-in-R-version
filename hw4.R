args = commandArgs(TRUE)
naivebayes = function(file1, file2, TAN) {
  options(warn = -1)
  options(digits = 12)
  #get dataset
  library(foreign)
  #train = read.arff("lymph_train.arff")
  #train = read.arff("vote_train.arff")
  train = read.arff(file1)
  #test = read.arff("lymph_test.arff")
  #test = read.arff("vote_test.arff")
  test = read.arff(file2)
  
  #get levels of each factor
  #getlevel = readLines("lymph_train.arff")
  #getlevel = readLines("vote_train.arff")
  getlevel = readLines(file1)
  gg = grep(x = getlevel, pattern = "attribute", value = T)
  alllevel = sub(x = gg, pattern = "@attribute \'.*\' \\{(.*)\\}", replacement = "\\1")
  alevel = gsub(x = alllevel, pattern = " |\'", replacement = "")
  levels = strsplit(x = alevel, split = ",")
  names(levels) = names(train)
  
  varnum = length(train)
  
  indepvar = levels[1:(varnum-1)]
  depvar = levels[varnum]
  
  
  ############################################
  #training part
  ############################################
  #get prob of y
  class0 = sum(train[,length(train)] == depvar[[1]][1])
  class1 = sum(train[,length(train)] == depvar[[1]][2])
  class = c(class0, class1)
  proby = c()
  for (classfac in class) {
    proby = c(proby, (classfac+1)/(dim(train)[1]+2))
  }
  names(proby) = c(depvar[[1]][1], depvar[[1]][2])
  
  
  #seperate train to different class level
  train1 = train[which(train[,length(train)] == depvar[[1]][1]),]
  train2 = train[which(train[,length(train)] == depvar[[1]][2]),]
  test1 = test[which(test[,length(test)] == depvar[[1]][1]),]
  test2 = test[which(test[,length(test)] == depvar[[1]][2]),]
  
  #without TAN
  
  #get frequencies of class0 and class1
  freq = function(train1, levels) {
    list1 = list()
    for (i in 1:length(train1)) {
      list1[[i]] = list()
    }
    names(list1) = names(train1)
    list1[length(train1)] = NULL
    for (i in 1:(varnum-1)) {
      levelnum = length(levels[[i]])
      for (j in 1:levelnum) {
        #judge = as.character(levels[[1]][[2]])
        judge = as.character(levels[[i]][[j]])
        #compare = as.character(train1[,1])
        compare = as.character(train1[,i])
        numhere = sum(compare == judge)
        probhere = round((numhere+1)/(dim(train1)[1]+levelnum), digits = 12)
        list1[[i]][[j]] = probhere
      }
      names(list1[[i]]) = names(levels[[i]])
    }
    return(list1)
  }
  
  class0ll1 = freq(train1, levels)
  class1ll1 = freq(train2, levels)
  
  #if naive bayes, print varnames and class
  if (TAN == "n") {
    for (varnamee in names(train)[1:(varnum-1)]) {
      cat(sep = " ", varnamee, "class", "\n")
    }
  }
  
  #with TAN
  if (TAN == "t") {
    #set a square matrix
    dataf = matrix(rep(0, (varnum-1)^2), nrow = (varnum-1))
    #get conditional mutual information of class0 and class1
    for (i in 1:length(indepvar)) {
      for (j in 1:length(indepvar)) {
        if (i != j) {
          for (k in 1:length(indepvar[[i]])) {
            for (o in 1:length(indepvar[[j]])) {
              lenvar1 = length(indepvar[[i]])
              lenvar2 = length(indepvar[[j]])
              nkk0 = sum(train1[,i] == indepvar[[i]][[k]] & train1[,j] == indepvar[[j]][[o]])
              nkk1 = sum(train2[,i] == indepvar[[i]][[k]] & train2[,j] == indepvar[[j]][[o]])
              pkk0 = (nkk0+1)/(dim(train)[1]+2*lenvar1*lenvar2)
              pkk1 = (nkk1+1)/(dim(train)[1]+2*lenvar1*lenvar2)
              n1 = dim(train1)[1]
              n2 = dim(train2)[1]
              pkky0 = (nkk0+1)/(dim(train1)[1]+lenvar1*lenvar2)
              pkky1 = (nkk1+1)/(dim(train2)[1]+lenvar1*lenvar2)
              pkyi0 = as.numeric(class0ll1[[i]][[k]])
              pkyj0 = as.numeric(class0ll1[[j]][[o]])
              pkyi1 = as.numeric(class1ll1[[i]][[k]])
              pkyj1 = as.numeric(class1ll1[[j]][[o]])
              ctb0 = round(pkk0 * (log(pkky0) - log(pkyi0) - log(pkyj0)), digits = 12) 
              ctb1 = round(pkk1 * (log(pkky1) - log(pkyi1) - log(pkyj1)), digits = 12)
              dataf[i,j] = dataf[i, j] + ctb0 + ctb1
            }
          }
        }
      }
    }
    
    #set Tree by using Prim's Algo
    dataparent = matrix(rep(FALSE, (varnum-1)^2), nrow = (varnum-1))
    passed = c(1)
    still = c(2:(varnum-1))
    while (length(still) > 0) {
      max = 0
      (maxpass = passed[1])
      (maxstill = still[1])
      for (s in still) {
        for (p in passed) {
          if (dataf[p, s] > max) {
            max = dataf[p, s]
            maxpass = p
            maxstill = s
          }
        }
      }
      dataparent[maxpass, maxstill] = TRUE
      passed = c(passed,maxstill)
      still = still[-which(still == maxstill)]
    }
    
    #print classes
    cat(sep = " ", names(train)[1], "class", "\n")
    for (i in 2: (varnum-1)) {
      parentcol = which(dataparent[,i] == TRUE)
      parent = names(train)[parentcol]
      here = names(train)[i]
      cat(sep = " ", here, parent, "class", "\n")
    }
      
    #update naive bayes
    for (i in 2:length(indepvar)) {
      parentnum = which(dataparent[, i] == TRUE)
      for (olo in 1:length(indepvar[[i]])) {
        for (klo in 1:length(indepvar[[parentnum]])) {
            nx11 = sum(train1[,parentnum] == indepvar[[parentnum]][[klo]] & train1[,i] == indepvar[[i]][[olo]])
            nx21 = sum(train1[,parentnum] == indepvar[[parentnum]][[klo]])
            nx12 = sum(train2[,parentnum] == indepvar[[parentnum]][[klo]] & train2[,i] == indepvar[[i]][[olo]])
            nx22 = sum(train2[,parentnum] == indepvar[[parentnum]][[klo]])
            lenelement = length(levels[[i]])
            lengthnow1 = length(class0ll1[[i]])
            numnum =  (nx11 + 1)/(nx21+lenelement)
            #self's level, parent's level, ratio given class = metastases
            class0ll1[[i]][[lengthnow1+1]] = c(olo, klo, numnum)
            
            lengthnow2 = length(class1ll1[[i]])
            numnum2 =  (nx12 + 1)/(nx22+lenelement)
            #self's level, parent's level, ratio given class = malign_lymph
            class1ll1[[i]][[lengthnow2+1]] = c(olo, klo, numnum2)
        }
      }
    }
  }
  
  ############################################
  #using test case to predict the data line by line
  ############################################
  countaccuracy = 0
  cat("\n")
  if (TAN == "t") {
    for (j in 1:dim(test)[1]) {
      #options(digits = 20)
      linehere = test[j,]
      #linehere = test[1,]
      roothere = names(table(linehere[1]))[which(table(linehere[1]) == 1)]
      levelindexroot = which(indepvar[[1]]==roothere)
      
      #print(linehere)
      #print(firsthere)
      #print(levelindexchild)
      #cat(sep = " ", linehere, firsthere, levelindexchild, "\n")
      
      #class = metastases
      rootphere0 = class0ll1[[1]][[levelindexroot]]
      classphere0 = proby[1]
      logprob0 = log(rootphere0 * classphere0)
      #class = malign_lymph
      rootphere1 = class1ll1[[1]][[levelindexroot]]
      classphere1 = proby[2]
      logprob1 = log(rootphere1 * classphere1)
      for (i in 2:(varnum-1)) {
        childtotallevel = sum(table(levels[[i]]))
        childhere = names(table(linehere[i]))[which(table(linehere[i]) == 1)]
        childindex = which(indepvar[[i]]==childhere)
        parent = which(dataparent[,i])
        parenttotallevel = sum(table(levels[[parent]]))
        parenthere = names(table(linehere[parent]))[which(table(linehere[parent]) == 1)]
        levelindexparent = which(indepvar[[parent]]==parenthere)
        logorigi0 = class0ll1[[i]][[childtotallevel + (childindex-1)*parenttotallevel + levelindexparent]][3]
        logorigi1 = class1ll1[[i]][[childtotallevel + (childindex-1)*parenttotallevel + levelindexparent]][3]
        logprob0 = logprob0 + log(logorigi0)
        logprob1 = logprob1 + log(logorigi1)  
      }
      
      prob0 = exp(logprob0)
      prob1 = exp(logprob1)
  
      if (prob0 > prob1) {
        #predict, real, prob
        acc = prob0/(prob0+prob1)
        #options(digits = 12)
        cat(sep = " ", levels[[varnum]][1], names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)], acc, "\n")
        if (levels[[varnum]][1] == names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)]) {
          countaccuracy = countaccuracy+1
        }  
      } else {
        acc = prob1/(prob0+prob1)
        options(digits = 12)
        cat(sep = " ", levels[[varnum]][2], names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)], acc, "\n")
        if (levels[[varnum]][2] == names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)]) {
          countaccuracy = countaccuracy+1
        }
      }
    }
  }
  if (TAN == "n") {
    for (j in 1:dim(test)[1]) {
      #options(digits = 18)
      linehere = test[j,]
      proclass0 = log(proby[1])
      proclass1 = log(proby[2])
      for (i in 1:(varnum-1)) {
        childhere = names(table(linehere[i]))[which(table(linehere[i]) == 1)]
        childindex = which(indepvar[[i]]==childhere)
        proclass0 = proclass0 + log(class0ll1[[i]][[childindex]])
        proclass1 = proclass1 + log(class1ll1[[i]][[childindex]])
      }
      prob0 = exp(proclass0)
      prob1 = exp(proclass1)
      if (prob0 > prob1) {
        #predict, real, prob
        acc = prob0/(prob0+prob1)
        options(digits = 12)
        cat(sep = " ", levels[[varnum]][1], names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)], acc, "\n")
        if (levels[[varnum]][1] == names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)]) {
          countaccuracy = countaccuracy+1
        }  
      } else {
        acc = prob1/(prob0+prob1)
        options(digits = 12)
        cat(sep = " ", levels[[varnum]][2], names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)], acc, "\n")
        if (levels[[varnum]][2] == names(table(linehere[varnum]))[which(table(linehere[varnum]) == 1)]) {
          countaccuracy = countaccuracy+1
        }
      }
    }
  }
  cat(sep = "", "\n", countaccuracy, "\n")
}

naivebayes(args[1], args[2], args[3])



