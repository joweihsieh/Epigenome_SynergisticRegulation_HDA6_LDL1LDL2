samtools merge SRR13443687.bam SRR13443687_1.bam SRR13443687_2.bam
samtools merge SRR13443688.bam SRR13443688_1.bam SRR13443688_2.bam &
samtools merge SRR24421255.bam SRR24421255_1.bam SRR24421255_2.bam &
samtools merge SRR24421259.bam SRR24421259_1.bam SRR24421259_2.bam &
samtools merge SRR24476325.bam SRR24476325_1.bam SRR24476325_2.bam &
samtools merge SRR24476326.bam SRR24476326_1.bam SRR24476326_2.bam &
samtools merge SRR5320832.bam SRR5320832_1.bam SRR5320832_2.bam &
samtools merge SRR5320833.bam SRR5320833_1.bam SRR5320833_2.bam &

#bsbolt Sort -I SRR13443687.bam -O SRR13443687_sorted.bam
bsbolt Sort -I SRR13443688.bam -O SRR13443688_sorted.bam &
bsbolt Sort -I SRR24421255.bam -O SRR24421255_sorted.bam &
bsbolt Sort -I SRR24421259.bam -O SRR24421259_sorted.bam &
bsbolt Sort -I SRR24476325.bam -O SRR24476325_sorted.bam &
bsbolt Sort -I SRR24476326.bam -O SRR24476326_sorted.bam &
bsbolt Sort -I SRR5320832.bam -O SRR5320832_sorted.bam &
bsbolt Sort -I SRR5320833.bam -O SRR5320833_sorted.bam &



bsbolt Sort -I SRR1569252_1.bam -O SRR1569252_sorted.bam &
bsbolt Sort -I SRR1569253_1.bam -O SRR1569253_sorted.bam & 
bsbolt Sort -I SRR5239972_1.bam -O SRR5239972_sorted.bam &
bsbolt Sort -I SRR5239973_1.bam -O SRR5239973_sorted.bam &
bsbolt Sort -I SRR5239984_1.bam -O SRR5239984_sorted.bam &
bsbolt Sort -I SRR5239985_1.bam -O SRR5239985_sorted.bam &
bsbolt Sort -I SRR5277784_1.bam -O SRR5277784_sorted.bam &
bsbolt Sort -I SRR5277785_1.bam -O SRR5277785_sorted.bam &
bsbolt Sort -I SRR6675996_1.bam -O SRR6675996_sorted.bam &
bsbolt Sort -I SRR6675997_1.bam -O SRR6675997_sorted.bam &
bsbolt Sort -I SRR7346214_1.bam -O SRR7346214_sorted.bam &
bsbolt Sort -I SRR7346216_1.bam -O SRR7346216_sorted.bam &

