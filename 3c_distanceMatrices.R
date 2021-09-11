library('Rcpp')

sourceCpp('include/bcdist_sparse.cpp')

## Run bray curtis distance on the sparse matrices.

args = commandArgs(trailingOnly = TRUE)
infile = args[1]
outfile = args[2]

if (tolower(stringi::stri_sub(outfile, from=-3)) != '.gz')
    outfile = paste0(outfile, '.gz')

## Read in SV table from infile.
svs = data.table::fread(infile, sep='\t',
    colClasses=c('character','character','integer'))
colnames(svs) = c('seqid', 'variable', 'value')

## pool svs unique to each sample into one row. This is a shortcut to reduce the
## immense number of calculations that we would otherwise need to perform.
nsams = svs[,.(N=.N), by='seqid']
uniquetons = nsams[N == 1, seqid]
uniquetonsvs = svs[seqid %in% uniquetons, .(seqid=paste0('singleton-', variable), value=sum(value)), by=variable][,.(seqid, variable, value)]
svpool = rbind(svs[!seqid %in% uniquetons,], uniquetonsvs)

## clear memory.
rm(uniquetons)
rm(svs)
rm(uniquetonsvs)

## Convert to integer matrix.
rows = svpool[,.(o=levels(as.factor(variable)))]
columns = svpool[,.(o=levels(as.factor(variable)))]
data.table::setkey(svpool, variable)

intcountdt = svpool[,.(x = as.numeric(as.factor(seqid)), y = as.numeric(as.factor(variable)), z = value)]
#spersam = intcountdt[,.(z=.N), by=y]

#spm = Matrix::sparseMatrix(i = intcountdt[,x], j = intcountdt[,y], x=intcountdt[,z])

tab = as.matrix(intcountdt)
#tab2 = as.matrix(spersam)

## Get combinations that will be performed.
sams = data.table::as.data.table(fastcombni(columns[,(1:.N)]))
bcd = data.table::data.table(a=columns[sams[,a],o],b=columns[sams[,b],o], d=bcdist(tab))
#bcd2 = data.table::data.table(a=columns[sams[,a],o],b=columns[sams[,b],o], d=bcdistsparse(tab, tab2))

data.table::fwrite(bcd, outfile, sep='\t')
