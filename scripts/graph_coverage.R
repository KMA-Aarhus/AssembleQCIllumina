suppressMessages(library("tidyverse"))

args <- commandArgs(TRUE)

path <- args[1]
#max_y <- as.integer(args[2])
#path <- "~/Documents/Projects/AssembleQC/A20_33780-KDD_N297-0335_S121_L001.coverage"

output_pth <- args[2] #str_replace(path, "\\.depth$","_coverage_mqc.jpg")
print(output_pth)
cvg <- read.table(path,header=FALSE)
colnames(cvg) = c("Contig","Pos","Reads")
large_contigs <- cvg %>% group_by(Contig) %>% count() %>% filter(n >= 1000) %>% pull(Contig)

filtered <- cvg %>% filter(Contig %in% large_contigs)

filtered <- cvg %>% mutate(contig_nr = as.numeric(factor(Contig))) %>% mutate(idx = row_number(), color_grp = contig_nr %% 9)
filtered <- filtered %>% mutate(color_grp = as.character(ifelse(color_grp==4,10, color_grp)))
reduced <- filtered %>%  slice(which(Pos %% 5 == 1))
below_thresh <- reduced %>% filter(Reads < 10 & Pos > 10)

plt <- ggplot(data=filtered, mapping=aes(x=idx, y=Reads, color=ifelse(Reads > 10, color_grp,"4")),show.legend=FALSE)
plt <- plt + geom_line(aes(group=Contig),show.legend=FALSE)
#plt <- plt + geom_text(data=below_thresh, aes(label=paste(Contig,Pos,sep=":")),color="Black",angle=90,check_overlap=TRUE,show.legend=FALSE)
plt <- plt + geom_hline(aes(yintercept = 10, linetype="10 Reads"))
plt <- plt + scale_color_brewer(palette="Paired")
plt <- plt + ylim(0,300) + labs(x="Index",y="# of Reads",linetype="Cutoff")

ggsave(output_pth,plot=plt, device="jpg")
