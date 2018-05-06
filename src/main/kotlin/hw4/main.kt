package hw4

import common.execCmd
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import kotlin.math.max
import kotlin.math.min

fun task1(outDir: String) {
    val dataDir = "/Johnny/students/NGS/data/4/E.coli"
    val testDataset = Dataset(
            "$dataDir/MG1655-K12.first10K.fasta",
            listOf("$dataDir/test_1.fastq",
                    "$dataDir/test_2.fastq")
    )

    val runDataset1 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            listOf("$dataDir/ecoli_mda_lane1_left.fastq.00.cor.fastq",
                    "$dataDir/ecoli_mda_lane1_right.fastq.00.cor.fastq")
    )

    val runDataset2 = Dataset(
            "$dataDir/MG1655-K12.fasta",
            listOf("$dataDir/s_6_1.fastq",
                    "$dataDir/s_6_2.fastq")
    )


    runTask1(testDataset, "test", outDir)

    val tasks = listOf(
            { runTask1(runDataset1, "e_coli_mda_lane", outDir) },
            { runTask1(runDataset2, "e_coli_s_6", outDir) }
    )
    tasks.parallelStream().forEach {
    //    it.invoke()
    }
}

fun task2(directory: String) {
    val dataDir = "/Johnny/students/NGS/data/5"
    val referenceName = "DH10B-K12.fasta"
    val referenceFile = "$directory/$referenceName"
    if (!File(referenceFile).exists()) {
        File("$dataDir/$referenceName").copyTo(File(referenceFile))
    }
    val testDataset = Dataset(
            referenceFile,
            listOf("$dataDir/tst.fastq")
    )

    val runDataset = Dataset(
            referenceFile,
            listOf("$dataDir/B22-730.fastq")
    )

    for (aligner in listOf(Aligner.BWA)) {
        runTask2(
                dataset = testDataset,
                runName = "DH10B_test_${aligner.name}",
                aligner = aligner,
                directory = directory
        )
        runTask2(
                dataset = runDataset,
                runName = "DH10B_run_${aligner.name}",
                aligner = aligner,
                directory = directory
        )
    }
}

fun main(args: Array<String>) {
    val outDir = "hw4_output"

    File(outDir).mkdirs()

    //task1(outDir)
    task2(outDir)
}