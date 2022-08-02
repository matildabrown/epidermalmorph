
nblue <- function(x){
  if(x==1) return("#0078D2")
  if(x==2) return(c("#004d88","#0090fe"))
  if(x==3) return(c("#004173","#0078D2","#69beff"))
  if(x==4) return(c("#009e91","#47fff0","#002088","#567eff"))
}

nred <- function(x){
  if(x==1) return("#F6511D")
  if(x==2) return(c("#ac2e07","#f76840"))
  if(x==3) return(c("#912706","#F6511D","#f98d6d"))
  if(x==4) return(c("#9f0748","#f74474","#ce4409","#f77239"))
}

nyellow <- function(x){
  if(x==1) return("#FFB400")
  if(x==2) return(c("#a67400","#ffbe26"))
  if(x==3) return(c("#8c6200","#FFB400","#ffd269"))
  if(x==4) return(c("#a65800","#ff9721","#a0a600","#ffb400"))
}

ngreen <- function(x){
  if(x==1) return("#7FB800")
  if(x==2) return(c("#527700","#45c964"))
  if(x==3) return(c("#456500","#7FB800","#b3f03a"))
  if(x==4) return(c("#657700","#c0e200","#2a7700","#37c83e"))
}

npurple <- function(x){
  if(x==1) return("#981867")
  if(x==2) return(c("#3f0a2b","#951765"))
  if(x==3) return(c("#350824","#981867","#d92293"))
  if(x==4) return(c("#5900ac","#b31dab","#3c0a18","#df3d6b"))
}
