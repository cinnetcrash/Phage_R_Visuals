# ==== Paketler ====
pkgs <- c("ape","genoPlotR","circlize")
inst <- pkgs[!(pkgs %in% rownames(installed.packages()))]
if(length(inst)) install.packages(inst, dependencies = TRUE)
library(ape); library(genoPlotR); library(circlize)

# ==== GFF yükle (FASTA'sız) ====
g1 <- read.gff("bakta_SA1298.nofasta.gff3")
g2 <- read.gff("bakta_SA1309.nofasta.gff3")

# ==== Sütunları normalize et ====
std_gff <- function(g){
  if ("type" %in% names(g) && !("feature" %in% names(g)))
    names(g)[names(g)=="type"] <- "feature"
  attr_col <- intersect(c("attribute","attributes","Attributes","ATTRIBUTES"), names(g))[1]
  if (is.na(attr_col)) g$attribute <- "" else names(g)[names(g)==attr_col] <- "attribute"
  need <- c("feature","start","end","strand","attribute")
  miss <- setdiff(need, names(g)); if(length(miss)) stop(paste("Eksik kolon:", paste(miss, collapse=", ")))
  g
}
g1 <- std_gff(g1); g2 <- std_gff(g2)

# ==== Kategorileme paleti ====
pal <- c(Structural="#2E86AB", Packaging="#E4572E", Replication="#76B041",
         Recomb="#A23B72", Lysis="#F1C40F", Hypo="#95A5A6", Other="#7B7FEC")

# ==== GFF -> dna_seg + renk (alt küme ile bire bir) ====
to_dna_colored <- function(g){
  keep <- subset(g, feature %in% c("gene","CDS","tRNA","rRNA"))
  # ürün metni
  prod <- ifelse(grepl("product=", keep$attribute, ignore.case=TRUE),
                 sub(".*product=([^;]+).*","\\1", keep$attribute, perl=TRUE), "")
  cat <- ifelse(grepl("capsid|portal|tail|collar|baseplate|fiber", prod, TRUE), "Structural",
                ifelse(grepl("terminase",  prod, TRUE), "Packaging",
                       ifelse(grepl("polymerase|primase|helicase|ligase", prod, TRUE), "Replication",
                              ifelse(grepl("integrase|recombinase|resolvase|transposase", prod, TRUE), "Recomb",
                                     ifelse(grepl("lysis|holin|endolysin", prod, TRUE), "Lysis",
                                            ifelse(grepl("hypothetical", prod, TRUE), "Hypo","Other"))))))
  colv <- pal[match(cat, names(pal))]
  
  # isim önceliği
  get_attr <- function(x, key){
    hit <- grepl(paste0(key,"="), x)
    out <- rep(NA_character_, length(x))
    out[hit] <- sub(paste0(".*",key,"=([^;]+).*"), "\\1", x[hit], perl=TRUE)
    out
  }
  nm_gene <- if ("gene" %in% names(keep)) keep$gene else NA_character_
  nm_ltag <- get_attr(keep$attribute, "locus_tag")
  nm_id   <- get_attr(keep$attribute, "ID")
  nm <- ifelse(!is.na(nm_gene)&nzchar(nm_gene), nm_gene,
               ifelse(!is.na(nm_ltag)&nzchar(nm_ltag), nm_ltag,
                      ifelse(!is.na(nm_id)&nzchar(nm_id), nm_id, keep$feature)))
  strand_num <- ifelse(keep$strand=="+", 1, -1)
  
  ds <- as.dna_seg(data.frame(name=nm, start=keep$start, end=keep$end,
                              strand=strand_num, stringsAsFactors=FALSE))
  ds$col <- colv              # uzunluklar bire bir
  list(ds=ds, keep=keep, cat=cat)
}

o1 <- to_dna_colored(g1); dna1 <- o1$ds
o2 <- to_dna_colored(g2); dna2 <- o2$ds

# ==== LINEER kontrol (isteğe bağlı) ====
plot_gene_map(
  dna_segs = list(dna1, dna2),
  dna_seg_labels = c("SA1298KAY-MASS9", "SA1309KAY-MASS15"),
  main = "Phage genomes (linear, colored by function)"
)

# ==== TEK DAİREDE İKİ GENOM + grid ====
normify <- function(ds){
  rng <- range(c(ds$start, ds$end)); len <- diff(rng)
  transform(as.data.frame(ds), x1=(start-rng[1])/len, x2=(end-rng[1])/len)
}
ds1 <- normify(dna1); ds2 <- normify(dna2)

circos.clear()
circos.par(track.height=0.08, gap.degree=2, start.degree=90, clock.wise=TRUE)
circos.initialize(factors="genome", xlim=c(0,1))

# 4 track: 1 outer +, 2 outer -, 3 inner +, 4 inner -
for(i in 1:4) circos.track("genome", ylim=c(0,1), track.height=0.08, bg.border=NA)

# blok çizici
draw_blocks <- function(ds, t_plus, t_minus){
  with(subset(ds, strand==1), circos.rect(xleft=x1, ybottom=0.15, xright=x2, ytop=0.95,
                                          col=col, border="black", track.index=t_plus))
  with(subset(ds, strand==-1), circos.rect(xleft=x1, ybottom=0.05, xright=x2, ytop=0.85,
                                           col=col, border="black", track.index=t_minus))
}
draw_blocks(ds1, 1, 2)  # SA1298 outer çift halka
draw_blocks(ds2, 3, 4)  # SA1309 inner çift halka

# Dış halka için grid ve cetvel (kb)
genome_len1 <- diff(range(c(dna1$start, dna1$end)))
step_kb <- max(2000, round(genome_len1/8/1000)*1000)  # ~8 ana tick
ticks <- seq(0, genome_len1, by=step_kb)
if (tail(ticks, 1) < genome_len1) ticks <- c(ticks, genome_len1)  # son tick'i ekle
major.at <- ticks / genome_len1

circos.trackPlotRegion("genome", track.index = 1, bg.border = NA,
                       panel.fun = function(x, y){
                         circos.axis(h="top", major.at = major.at,
                                     labels = paste0(ticks/1000, " kb"),
                                     labels.cex = 0.6, major.tick.length = convert_y(2, "mm"))
                       })

# Genom etiketleri
circos.text(-0.5, 1.2, "SA1298KAY-MASS9 (outer ±)", sector.index="genome", track.index=2, cex=0.9, facing="bending.outside")
circos.text(-0.5, 1.2, "SA1306KAY-MASS15 (inner ±)", sector.index="genome", track.index=4, cex=0.9, facing="bending.outside")
# Lejand
legend("topright",
       fill = pal, legend = names(pal), bty="n", cex=0.8, ncol=1,
       title = "Functional categories")

# ==== (İsteğe bağlı) GC içerik / GC skew trackleri ====
# FASTA yoksa atlar. Aynı klasörde SA1298.fasta ve SA1309.fasta varsa etkinleşir.
add_gc_tracks <- function(fasta_path, ring_index, window = 500){
  if (!requireNamespace("seqinr", quietly = TRUE)) return(invisible(NULL))
  seq <- seqinr::read.fasta(fasta_path, as.string=TRUE)[[1]]
  s <- unlist(strsplit(toupper(seq), ""))
  n <- length(s)
  w <- window
  gc <- gcskew <- rep(NA_real_, n-w+1)
  for(i in 1:(n-w+1)){
    win <- s[i:(i+w-1)]
    g <- sum(win=="G"); c <- sum(win=="C"); a <- sum(win=="A"); t <- sum(win=="T")
    gc[i] <- (g+c)/w
    gcskew[i] <- ifelse((g+c)==0, 0, (g-c)/(g+c))
  }
  x <- seq(0, 1, length.out=length(gc))
  circos.track("genome", ylim=c(min(gcskew,na.rm=TRUE), max(gcskew,na.rm=TRUE)),
               track.height=0.06, bg.border=NA)
  circos.lines(x, gcskew, sector.index="genome", track.index=ring_index, col="#333333")
}

# Örnek kullanım (dosyalarınız varsa açın):
add_gc_tracks("SA1298KAY-MAS-S9.fasta", ring_index = 5)
add_gc_tracks("SA1306-KAY-MAS-S15.fasta", ring_index = 6)

# Kayıt için:
png("both_genomes_single_circle.png", 1800, 1800, res=300);  ## sonra tüm çizmeleri tekrar çalıştırın
dev.off()
