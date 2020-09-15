
################################################################################
#    &&&....&&&    % Project: Mapping KO ids to gene symbols  #
#  &&&&&&..&&&&&&  % Authors: Xiner Nie, Bo Li, Jingxin Tao, Hao He            #
#  &&&&&&&&&&&&&&  % Date: Aug. 18th, 2020                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.6.0;                             #
#       &&&&       % Platform: x86_64-w64-mingw32/x64 (64-bit)                 #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 10: Evaluation Criteria for data clustering.
### ****************************************************************************

download.file("http://rest.kegg.jp/list/ko", "ko.txt")

anno <- readLines("ko.txt")

strsplit(anno[1], "\t")[[1]][1]
strsplit(anno[1], "\t")[[1]][2]

strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]]

strsplit(strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]][2], " \\[")[[1]][1]

gsub("EC", "\\[EC", strsplit(strsplit(strsplit(anno[1], "\t")[[1]][2], "; ")[[1]][2], " \\[")[[1]][2])












