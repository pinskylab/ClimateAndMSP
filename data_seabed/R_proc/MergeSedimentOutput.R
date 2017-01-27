### Compile sediment data into single datatable
setwd("~/Documents/Collaborations/Rutgers/Sediment Grain")

library(data.table)

gmx <- as.data.table(readRDS("haul_sp/hauls_sed_gmx_db.rds"))
pac <- as.data.table(readRDS("haul_sp/hauls_sed_pac_idw.rds"))
ak <- as.data.table(readRDS("haul_sp/hauls.ak.idw.rds"))
atl <- as.data.table(readRDS("haul_sp/hauls.atl.sed.conmap.rds"))
lab <- as.data.table(readRDS("haul_sp/hauls.lab.idw.rds"))
scot <- as.data.table(readRDS("haul_sp/hauls.scot.idw.rds"))
newf <-as.data.table(readRDS("haul_sp/hauls.newf.poly.rds"))
sogulf <- as.data.table(readRDS("haul_sp/hauls.sogulf.poly.rds"))

atl[,"SEDIMENT":=NULL]
atl[,"SEDNUM":=NULL]
lab[,"lon.adj":=lon]
scot[,"lon.adj":=lon]
newf[,"lon.adj":=lon]

### Add info on which method was used to get haul sediment information
gmx[,"method":="poly.extract"]
atl[,"method":="poly.extract"]
newf[,"method":="poly.extract"]
sogulf[,"method":="poly.extract"]
ak[,"method":="idw"]
pac[,"method":="idw"]
scot[,"method":="idw"]
lab[,"method":="idw"]



sediment.all <- rbind(gmx, pac, ak, atl, lab, scot, newf, sogulf)
saveRDS(sediment.all, "haul_sp/sediment.all.out.rds")

### hauls that were interpolated twice-- only on ScotianShelfSpring 
## because hauls overlapped with two sediment datasets (once with polygon extract of conmap and second with idw of Expedition database)
sediment.dup <- sediment.all[duplicated(haulid)]
sediment.all[haulid=="HAM1979013-005"]