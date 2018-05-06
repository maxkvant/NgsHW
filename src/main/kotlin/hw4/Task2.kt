package hw4

import java.io.File

fun runTask2(dataset: Dataset, runName: String, aligner: Aligner, directory: String) {
    val samFile = "$directory/${runName}_aln.sam"

    if (!File(samFile).exists()) {
        when (aligner) {
            Aligner.Bowtie2 -> runBowtie2(dataset, samFile)
            Aligner.BWA -> runBwa(dataset, samFile)
        }
    }
    val processor = DatasetProcessor(
            dataset = dataset,
            runName = runName,
            samFile = samFile,
            directory = directory
    )

    processor.plotCoverage()
    processor.countSubstitutionTable()
    processor.plotMonomers()
}