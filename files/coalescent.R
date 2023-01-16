# Draw a tree

mrca.id = function(n) {
  2*n-1
}

new.tree = function(n) {
  cols = c("id","dec1","dec2","anc","time","length")
  mrca = mrca.id(n)
  t = matrix(0,mrca,length(cols))
  colnames(t) = cols
  t[,"id"] = 1:mrca
  return(t)
}

sim.coalescent.tree = function(n, PNe=1) {
  tree = new.tree(n)
  active = 1:n
  time = 0
  tree[active,"time"] = time
  mrca = mrca.id(n)
  for(i in (n+1):mrca) {
    k = mrca-i+2
    aid = sample(1:k,2)
    coal = active[aid]
    dt = rexp(1,choose(k,2))
    time = time + dt
    tree[coal,"anc"] = i
    tree[i,c("dec1","dec2")] = coal
    tree[i,"time"] = time
    tree[coal,"length"] = time-tree[coal,"time"]
    keep = rep(TRUE,k)
    keep[aid[2]] = FALSE
    active[aid[1]] = i
    active = active[keep]
  }
  tree[,c("time","length")] = tree[,c("time","length")]*PNe
  return(tree)
}

count.descendants = function(tree) {
  mrca = nrow(tree)
  n = (mrca+1)/2
  ndec = rep(0,mrca)
  for(i in (n+1):mrca) {
    dec0 = tree[i,"dec1"]
    dec1 = tree[i,"dec2"]
    ndeci = ndec[c(dec0,dec1)]
    ndeci[ndeci==0] = 1
    ndec[i] = sum(ndeci)
  }
  return(ndec)
}

draw.tree.skewwhiff = function(tree,ylim=range(tree[,"time"])) {
  mrca = nrow(tree)
  n = (mrca+1)/2
  plot(c(1,n),ylim,type="n",axes=F,ann=F)
  axis(2)
  tree = as.matrix(cbind(as.data.frame(tree),"ndec"=count.descendants(tree),"x"=n/2))
  tree[tree[,"ndec"]==0,"ndec"] = 1
  for(i in mrca:(n+1)) {
    xpoint = tree[i,"x"]
    xwidth = tree[i,"ndec"]
    xlo = xpoint-xwidth/2
    xhi = xpoint+xwidth/2
    ypoint = tree[i,"time"]
    dec = tree[i,c("dec1","dec2")]
    xdec = tree[dec,"ndec"]
    xdec = xlo + cumsum(xdec)-xdec+xdec/2
    ydec = tree[dec,"time"]
    lines(xdec,rep(ypoint,2))
    lines(rep(xdec[1],2),c(ydec[1],ypoint))
    lines(rep(xdec[2],2),c(ydec[2],ypoint))
    tree[dec,"x"] = xdec
  }
}


draw.tree.simple = function(tree,ylim=range(tree[,"time"]),padding=1,balance="adhoc",...) {
  balance = pmatch(balance,c("adhoc","left","right","sample"))
  if(is.na(balance)) stop("balance must be one of 'adhoc' 'left' 'right' or 'sample'")
  mrca = nrow(tree)
  n = (mrca+1)/2
  plot(c(1-padding/2,n+padding/2),ylim,type="n",axes=F,ann=F)
  axis(2)
  tree = as.matrix(cbind(as.data.frame(tree),"ndec"=count.descendants(tree),"x"=n/2,"x2"=0))
  tree[tree[,"ndec"]==0,"ndec"] = 1
  for(i in mrca:(n+1)) {
    xpoint = tree[i,"x"]
    xwidth = tree[i,"ndec"]
    xlo = xpoint-xwidth/2
    xhi = xpoint+xwidth/2
    ypoint = tree[i,"time"]
    dec = tree[i,c("dec1","dec2")]
    xdec = tree[dec,"ndec"]
    if(balance!=1) {
      if(balance==4) dec=sample(dec)
      else dec=dec[order(xdec,decreasing=(balance==2))]
    }
    xdec = tree[dec,"ndec"]
    xdec = xlo + cumsum(xdec)-xdec+xdec/2
    ydec = tree[dec,"time"]
    tree[dec,"x"] = xdec
  }
  tree[1:n,"x2"] = rank(tree[1:n,"x"])
  for(i in 1:n) lines(rep(tree[i,"x2"],2),tree[i,"time"]+c(0,tree[i,"length"]),...)
  for(i in (n+1):(mrca-1)) {
    ypoint = tree[i,"time"]
    dec = tree[i,c("dec1","dec2")]
    xdec = tree[dec,"x2"]
    lines(xdec,rep(ypoint,2),...)
    xbar = mean(xdec)
    lines(rep(xbar,2),tree[i,"time"]+c(0,tree[i,"length"]),...)
    tree[i,"x2"] = xbar
  }
  ypoint = tree[mrca,"time"]
  dec = tree[mrca,c("dec1","dec2")]
  xdec = tree[dec,"x2"]
  lines(xdec,rep(ypoint,2),...)
}

sim.serial.coalescent.tree = function(n, PNe=1, times=rep(0,n),expectedtimes=FALSE) {
  if(length(PNe)!=1) stop("PNe must be a scalar")
  if(length(times)==1) times = rep(times,n)
  if(length(times)!=n) stop("times must have length 1 or n")
  times = sort(times)
  tree = new.tree(n)
  time = min(times)
  nextsample = suppressWarnings(min(times[times>time]))
  active = which(times==time)
  tree[1:n,"time"] = times
  mrca = mrca.id(n)
#  while((k=length(active))>1 | nextsample<Inf) {
  for(i in (n+1):mrca) {
    rt = ifelse(expectedtimes,0.5,runif(1))
    while(TRUE) {
      k = length(active)
	  rate = choose(k,2)
#      dt = ifelse(k>1,ifelse(expectedtimes,1/choose(k,2),rexp(1,choose(k,2))),Inf)
      dt = ifelse(k>1,qexp(rt,rate),Inf)
      if(time+dt>=nextsample) {
	    rt = ifelse(expectedtimes,rt - pexp(nextsample-time,rate),runif(1))
        time = nextsample
        active = c(active,which(times==time))
        nextsample = suppressWarnings(min(times[times>time]))
      }
      else break
    }
    aid = sample(1:k,2)
    coal = active[aid]
    time = time + dt
    tree[coal,"anc"] = i
    tree[i,c("dec1","dec2")] = coal
    tree[i,"time"] = time
    tree[coal,"length"] = time-tree[coal,"time"]
    keep = rep(TRUE,k)
    keep[aid[2]] = FALSE
    active[aid[1]] = i
    active = active[keep]
  }
  tree[,c("time","length")] = tree[,c("time","length")]*PNe
  return(tree)
}

iam.mutate.tree = function(tree,theta,expected=FALSE) {
  totalrate = sum(tree[,"length"])*theta/2
  nmut = ifelse(expected,round(totalrate),rpois(1,totalrate))
  cols = c("id","state","branch","time")
  muts = matrix(0,nrow=nmut+1,ncol=length(cols))
  colnames(muts) = cols
  branch = c(rep(1:nrow(tree),times=rmultinom(1,nmut,tree[,"length"])),nrow(tree))
  muts[,"branch"] = branch
  muts[,"time"] = tree[branch,"time"]+tree[branch,"length"]*runif(nmut+1)
  ord = order(muts[,"branch"],muts[,"time"])
  muts[ord,"id"] = 1:(nmut+1)
  muts[ord,"state"] = (nmut+1):1
  return(muts[ord,,drop=FALSE])
}

draw.tree = function(tree,ylim=range(tree[,"time"]),padding=1,balance="adhoc",plot.it=TRUE,...) {
  balance = pmatch(balance,c("adhoc","left","right","sample"))
  if(is.na(balance)) stop("balance must be one of 'adhoc' 'left' 'right' or 'sample'")
  mrca = nrow(tree)
  n = (mrca+1)/2

  tree = as.matrix(cbind(as.data.frame(tree),"ndec"=count.descendants(tree),"x"=n/2,"x2"=0))
  tree[tree[,"ndec"]==0,"ndec"] = 1
  for(i in mrca:(n+1)) {
    xpoint = tree[i,"x"]
    xwidth = tree[i,"ndec"]
    xlo = xpoint-xwidth/2
    xhi = xpoint+xwidth/2
    ypoint = tree[i,"time"]
    dec = tree[i,c("dec1","dec2")]
    xdec = tree[dec,"ndec"]
    if(balance!=1) {
      if(balance==4) dec=sample(dec)
      else dec=dec[order(xdec,decreasing=(balance==2))]
    }
    xdec = tree[dec,"ndec"]
    xdec = xlo + cumsum(xdec)-xdec+xdec/2
    ydec = tree[dec,"time"]
    tree[dec,"x"] = xdec
  }
  tree[1:n,"x2"] = rank(tree[1:n,"x"])

  cols = c("id","node","x1","x2","y1","y2")
  LINES = matrix(0,nrow=3*n-3,ncol=length(cols))
  colnames(LINES) = cols
  ctr = 1
  for(i in 1:n) {
    LINES[ctr,] = c(ctr,i,rep(tree[i,"x2"],2),tree[i,"time"]+c(0,tree[i,"length"]))
	ctr = ctr+1
  }
  for(i in (n+1):(mrca-1)) {
    ypoint = tree[i,"time"]
    dec = tree[i,c("dec1","dec2")]
    xdec = tree[dec,"x2"]
    LINES[ctr,] = c(ctr,i,xdec,rep(ypoint,2))
    xbar = mean(xdec)
	LINES[ctr+1,] = c(ctr+1,i,rep(xbar,2),tree[i,"time"]+c(0,tree[i,"length"]))
    tree[i,"x2"] = xbar
	ctr = ctr+2
  }
  ypoint = tree[mrca,"time"]
  dec = tree[mrca,c("dec1","dec2")]
  xdec = tree[dec,"x2"]
  LINES[ctr,] = c(ctr,mrca,xdec,rep(ypoint,2))
  if(plot.it) {
    plot(c(1-padding/2,n+padding/2),ylim,type="n",axes=F,ann=F)
    axis(2)
    for(i in 1:nrow(LINES)) lines(LINES[i,c("x1","x2")],LINES[i,c("y1","y2")],...)
  }
  else return(LINES)
}

draw.lines = function(LINES,...) {
  for(i in 1:nrow(LINES)) lines(LINES[i,c("x1","x2")],LINES[i,c("y1","y2")],...)
}

draw.mutated.tree = function(tree,mutations,ylim=range(tree[,"time"]),padding=1,balance="adhoc",plot.it=TRUE,clr=NULL,...) {
  LINES = draw.tree(tree,ylim,padding,balance,plot.it=FALSE)
  
  # Need to split vertical nodes according to number of mutations
  # Each mutation adds a vertical node except at the root
  cols = c(colnames(LINES),"state")
#  nrow(LINES) = n+  2*(n-1)-1 = 3*n-3
  n = (nrow(LINES)+3)/3
  nmut = nrow(mutations)-1
  nr = nrow(LINES)+nrow(mutations)-1
  CLINES = matrix(0,nrow=nr,ncol=length(cols))
  colnames(CLINES) = cols
  
  # Choose colour scheme
  if(is.null(clr)) {
    if(min(mutations[,"state"])<1) stop("Mutations must be numbered positively")
    clr = sample(rainbow(max(mutations[,"state"])))
  }
  RGB = matrix(0,nrow=length(clr),ncol=3)
  colnames(RGB) = c("r","g","b")
  for(i in 1:length(clr)) {
    RGB[i,] = as.numeric(paste("0x",c(substr(clr[i],2,3),substr(clr[i],4,5),substr(clr[i],6,7)),sep=""))
  }
  
  # il indexes LINES, ic indexes CLINES
  # im indexes mutations and it indexes the tree
  il = nrow(LINES); ic = nr; im = nrow(mutations)
  CLINES[ic,] = c(LINES[il,],mutations[im,"state"])
  il = il-1; ic = ic-1; im = im-1
  # Go through the non-root coalescent nodes
  for(it in (2*n-2):1) {
    ypoint = LINES[il,"y2"]
	anc = tree[it,"anc"]
	canc = ic-1+which(CLINES[ic:nr,"node"]==anc)[1]
	amut = CLINES[canc,"state"]
	while(TRUE) {
	  if(im==0) break
	  if(mutations[im,"branch"]!=it) break
	  CLINES[ic,] = c(LINES[il,c("id","node","x1","x2")],mutations[im,"time"],ypoint,amut)
	  ypoint = mutations[im,"time"]
	  amut = mutations[im,"state"]
	  ic = ic-1; im = im-1
	}
	CLINES[ic,] = c(LINES[il,c("id","node","x1","x2","y1")],ypoint,amut)
	if(it>n) {
      CLINES[ic-1,] = c(LINES[il-1,],amut)
      ic = ic-2; il = il-2
	}
	else {
      ic = ic-1; il = il-1
	}
  }
  CLINES[,"id"] = 1:nr
  
  if(plot.it) {
    plot(c(1-padding/2,n+padding/2),ylim,type="n",axes=F,ann=F)
    axis(2)
    for(i in 1:nrow(CLINES)) lines(CLINES[i,c("x1","x2")],CLINES[i,c("y1","y2")],col=clr[CLINES[i,"state"]],...)
  }
  else return(CLINES)
}

superimpose.mutations = function(tree,mutations,...) {
  draw.tree(tree)
  
}

genotype = function(tree,mutation) {
  state = rep(0,nrow(tree))
  state[nrow(tree)] = 0
  wh = which(mutation[,"branch"] == nrow(tree))
  if(length(wh)>0) state[nrow(tree)] = max(mutation[wh,"state"])
  for(i in (length(state)-1):1) {
    state[i] = state[tree[i,"anc"]]
	wh = which(mutation[,"branch"] == i)
	if(length(wh)>0) state[i] = max(mutation[wh,"state"])
  }
  return(state)
}

# Convert tree times to exponential growth (r is positive for a population shrinking into the past)
expgrowth = function(tree,r) {
	tree[,"time"] = 1/r*log(1+r*tree[,"time"])
	ix = 1:(nrow(tree)-1)
	tree[ix,"length"] = tree[tree[ix,"anc"],"time"]-tree[ix,"time"]
	return(tree)
}
#while(TRUE) {draw.tree(expgrowth(sim.coalescent.tree(49),-2)); Sys.sleep(.5);}

# Calculate a doubling time from a growth rate
dtime = function(r) 1/r*log(2)
