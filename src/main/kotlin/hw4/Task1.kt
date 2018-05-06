package hw4

import java.io.File

fun runTask1(dataset: Dataset, runName: String, directory: String) {
    val samFile = "$directory/${runName}_aln.sam"
    if (!File(samFile).exists()) {
        runBowtie2(dataset, samFile)
    }
    val processor = DatasetProcessor(
            dataset = dataset,
            runName = runName,
            samFile = samFile,
            directory = directory
    )

    processor.plotCoverage()

    val tasks = listOf(
            { processor.countAlignedReads() },
            { processor.plotInsertSizes() },
            { processor.countSubstitutionTable() }
    )
    tasks.parallelStream().forEach {
        it.invoke()
    }
}