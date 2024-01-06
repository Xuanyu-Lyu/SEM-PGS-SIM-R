
#Automates Start Values for Mx

st.parameters <- c("a","b","d","e","s","t","f","vt")
st.values <- round(runif(8,.1,.9),3)   
st.elements <- c(" D 1 1 E 1 1",
                 " F 1 1",
                 " N 1 1 O 1 1",
                 " W 1 1 X 1 1",
                 " Q 1 1 R 1 1",
                 " T 1 1 U 1 1",
                 " J 1 1 K 1 1",
                 " A 1 1 B 1 1 C 1 1 D 1 1")
     
st.mat <- matrix(paste("#define $st_",st.parameters, "            Start ",st.values, st.elements,sep=""),nrow=length(st.values))

write.table(st.mat,file="st.val.mat",append=FALSE,row.names=FALSE,col.names=FALSE,quote=FALSE)

