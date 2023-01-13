### Read a FASTA file
read.fasta <- function(fname, as.char=FALSE) {
  a = scan(fname,what=character(0),sep="\n",quiet=TRUE,na.strings="")
  wh = as.vector(sapply(a,substr,1,1))==">"
  labs = substr(as.character(a[wh]),2,1000);

  lseqs = a[!wh]
  nlines = length(lseqs)%/%length(labs)
  n = length(lseqs)%/%nlines
  seqs = rep("",n);
  names(seqs) <- labs
  for(i in 1:n) {
    ibeg = (i-1)*nlines+1
    iend = i*nlines
    seqs[i] = paste(lseqs[ibeg:iend],collapse="")
  }
  seqlen = as.numeric(sapply(seqs,nchar))
  if(length(seqlen)>1 & var(seqlen)>0) {
    warning("Sequences have differing lengths");
	mx = max(seqlen)
	for(i in 1:n) seqs[i] = paste(seqs[i],paste(rep("-",mx-seqlen[i]),collapse=""),sep="")
  }
  L = as.numeric(nchar(seqs[1]))

  SEQ = array("-",dim=c(n,L))
  for(i in 1:n) SEQ[i,] = unlist(strsplit(seqs[i],""))
  rownames(SEQ) <- labs;
  if(as.char==TRUE) {
    return(SEQ);
  } else {
    fSEQ = apply(toupper(SEQ),2,factor,levels=c("A","G","C","T"));
    return(fSEQ);
  }
}

### Write a FASTA file
write.fasta <- function(DNA,filename) {
  ofile <- file(filename,"w");
  for(n in 1:nrow(DNA)) {
    writeLines(paste(">",rownames(DNA)[n],sep=""),ofile);
    writeLines(paste(DNA[n,],collapse=""),ofile);    
  }
  close(ofile);
}

### General
totriplet = function(x) {
  L = floor(length(x)/3)*3
  paste(x[seq(1,L,by=3)],x[seq(2,L,by=3)],x[seq(3,L,by=3)],sep="")
}
geneticCode = list(
"TTT"="Phe","TTC"="Phe","TTA"="Leu","TTG"="Leu",
"TCT"="Ser","TCC"="Ser","TCA"="Ser","TCG"="Ser",
"TAT"="Tyr","TAC"="Tyr","TAA"="STO","TAG"="STO",
"TGT"="Cys","TGC"="Cys","TGA"="STO","TGG"="Trp",
"CTT"="Leu","CTC"="Leu","CTA"="Leu","CTG"="Leu",
"CCT"="Pro","CCC"="Pro","CCA"="Pro","CCG"="Pro",
"CAT"="His","CAC"="His","CAA"="Gln","CAG"="Gln",
"CGT"="Arg","CGC"="Arg","CGA"="Arg","CGG"="Arg",
"ATT"="Ile","ATC"="Ile","ATA"="Ile","ATG"="Met",
"ACT"="Thr","ACC"="Thr","ACA"="Thr","ACG"="Thr",
"AAT"="Asn","AAC"="Asn","AAA"="Lys","AAG"="Lys",
"AGT"="Ser","AGC"="Ser","AGA"="Arg","AGG"="Arg",
"GTT"="Val","GTC"="Val","GTA"="Val","GTG"="Val",
"GCT"="Ala","GCC"="Ala","GCA"="Ala","GCG"="Ala",
"GAT"="Asp","GAC"="Asp","GAA"="Glu","GAG"="Glu",
"GGT"="Gly","GGC"="Gly","GGA"="Gly","GGG"="Gly")
oneLetterCodes = unlist(list("Ala"="A","Arg"="R","Asn"="N","Asp"="D","Cys"="C","Glu"="E","Gln"="Q","Gly"="G","His"="H","Ile"="I","Leu"="L","Lys"="K","Met"="M","Phe"="F","Pro"="P","Ser"="S","Thr"="T","Trp"="W","Tyr"="Y","Val"="V","STO"="X","---"="-"))
aminoAcids = names(table(unlist(geneticCode)))
oneLetterAminoAcids = names(table(unlist(oneLetterCodes)))
tripletNames = names(geneticCode)

transcribe = function(x) {
  y = t(sapply(1:nrow(x),function(i) totriplet(x[i,])))
  rownames(y) = rownames(x)
  return(y)
}
translate = function(x,oneLetter=FALSE) {
  x = toupper(x)
  tr = t(apply(x,1,function(y)sapply(y,function(i) {aa=geneticCode[[i]];ifelse(is.null(aa),"---",aa)} )))
  if(oneLetter) tr = t(apply(tr,1,function(y) oneLetterCodes[y]))
  rownames(tr) = rownames(x)
  return(tr)
}
view.nucleotide = function(x) {
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x,levels=c("-","A","G","C","T"))),nrow=nrow(x))),col=c("white","red","green","yellow","blue"))
}
view.codon = function(x) {
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x),levels=tripletNames),nrow=nrow(x))),col=rainbow(20))
}
view.protein = function(x,oneLetter=FALSE) {
  levs = aminoAcids
  if(oneLetter) levs = oneLetterAminoAcids
  cols = rainbow(20)
  if(oneLetter) cols = c("white",cols)
  image(0:ncol(x),0:nrow(x),t(matrix(as.numeric(factor(x,levels=levs)),nrow=nrow(x))),col=cols)
}
# Assumes a fasta file representing a single genome, possibly split across contigs
read.fasta.ref = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	beg = substr(r,1,1)
	gd = beg!=">"
	rcat = paste(r[gd],collapse="")
	return(toupper(unlist(strsplit(rcat,""))))
}
# Assumes a fasta file representing a single genome, possibly split across contigs
read.fasta.ref.contig = function(ref_file) {
	r = scan(ref_file,what=character(0),sep="\n")
	beg = substr(r,1,1)
	gd = beg!=">"
	contig = rep(cumsum(!gd)[gd],times=nchar(r[gd]))
	return(contig)
}

