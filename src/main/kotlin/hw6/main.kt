package hw6

import common.execCmd
import common.forSamRecordsSingle
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequence
import java.io.File

const val starPath = "/Johnny/students/NGS/data/old/7/STAR-2.6.0b/bin/Linux_x86_64/STAR"

const val spadesPath = "/home/m_vinnichenko/intership/algorithmic-biology/assembler/spades.py"

const val rnaQuastPath = "/Johnny/students/NGS/data/8/rnaQUAST-1.5.1/rnaQUAST.py"

const val numThreads = 8

class Dataset(val reference: String, val reads: List<String>, val genesFile: String)


fun runStar(dataset: Dataset, outDir: String): List<String> {
    File(outDir).mkdirs()
    fun getSamFiles() = File(outDir).listFiles()
            .filter { it.isFile && it.name.endsWith(".sam") }
            .map { it.absolutePath }

    if (getSamFiles().isEmpty()) {
        val indexDir = File(outDir).absolutePath


        val cmdIndex = "$starPath --runThreadN $numThreads " +
                "--runMode genomeGenerate " +
                "--genomeDir ${indexDir} " +
                "--genomeFastaFiles ${dataset.reference} " +
                "--sjdbGTFfile ${dataset.genesFile}"
        execCmd(cmdIndex)

        val cmdRun = "$starPath --runThreadN $numThreads " +
                "--genomeDir $indexDir/ " +
                "--readFilesIn ${dataset.reads.joinToString(" ")} " +
                "--outFileNamePrefix $outDir/star"
        execCmd(cmdRun)
    }
    return getSamFiles()
}

fun task2(dataset: Dataset, outDir: String) {
    val samFiles = runStar(dataset, outDir)
    require(samFiles.size == 1)
    val samFile = samFiles[0]

    data class Segment(val l: Int, val r: Int, val referenceName: String)

    val geneExons = mutableMapOf<String, MutableList<Segment>>()

    File(dataset.genesFile).forEachLine { line ->
        if (line.startsWith("#") || line.isEmpty()) {
            return@forEachLine
        }
        val chunks = line.split("\t")
        val sequence = chunks[0]
        val type = chunks[2]
        val start = chunks[3].toInt()
        val end = chunks[4].toInt()
        val comment = chunks[8].split(";")[0]
        if (comment.startsWith("gene_id") && type == "exon") {
            if (!geneExons.containsKey(comment)) {
                geneExons[comment] = mutableListOf()
            }
            val segment = if (start <= end) Segment(start, end, sequence) else { print("!"); Segment(end, start, sequence) }
            geneExons[comment]!!.add(segment)
        }
    }

    val sequenceDict: Map<String, String> = FastaSequenceFile(File(dataset.reference), false).use { fastaFile ->
        val sequences = mutableListOf<ReferenceSequence>()
        var sequence: ReferenceSequence? = fastaFile.nextSequence()
        while (sequence != null) {
            sequences.add(sequence)
            sequence = fastaFile.nextSequence()
        }
        sequences.map {
            Pair(it.name.split(" ")[0], it.baseString)
        }.toMap()
    }

    val coverage: Map<String, IntArray> = sequenceDict
            .mapValues { IntArray(it.value.length) }
            .toMap()


    forSamRecordsSingle(samFile, false) { samRecord ->
        for (block in samRecord.alignmentBlocks) {
            for (pos in 0 until block.length) {
                val posReference = block.referenceStart - 1 + pos
                coverage[samRecord.referenceName]!![posReference] += 1
            }
        }
    }

    val minCoverage = 3
    var coveredGenes = 0

    for (pGeneExons in geneExons) {
        val positions = mutableSetOf<Int>()
        val positionsCovered = mutableSetOf<Int>()

        for (segment in pGeneExons.value) {
            for (pos in segment.l..segment.r) {
                if (coverage[segment.referenceName]!![pos] >= minCoverage) {
                    positionsCovered.add(pos)
                }
                positions.add(pos)
            }
        }
        if (positionsCovered.size >= positions.size * 0.95) {
            coveredGenes += 1
        }
    }

    val genesCnt = geneExons.size
    println(">= 95% coverage: $coveredGenes of $genesCnt (${coveredGenes / genesCnt.toDouble() * 100}%)")

}

fun task3(dataset: Dataset, outDir: String) {
    val assembleDir = "$outDir/assemble"
    if (!File(assembleDir).exists()) {
        File(assembleDir).mkdirs()

        require(dataset.reads.size == 2)

        val spadesCmd = "$spadesPath --rna " +
                "-1 ${dataset.reads[0]} -2 ${dataset.reads[1]} " +
                "--threads $numThreads " +
                "-o $outDir "
        execCmd(spadesCmd)
    }

    val quastCmd = "python2 $rnaQuastPath --transcripts $assembleDir/transcripts.fasta " +
            "-r ${dataset.reference} " +
            "--gtf ${dataset.genesFile} " +
            "-o $outDir/quast "
    execCmd(quastCmd)
}

fun main(args: Array<String>) {
    val outDir = "hw6_output"
    File(outDir).mkdirs()
    val dataDir = "/Johnny/students/NGS/data/8"


    val dataset = Dataset(
            reference = "$dataDir/ref.fa",
            reads = listOf("$dataDir/SRR453566_1.fastq", "$dataDir/SRR453566_2.fastq"),
            genesFile = "$dataDir/genes.gtf"
    )

    task2(dataset, outDir)
    task3(dataset, outDir)
}