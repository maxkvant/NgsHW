package hw1

import com.github.sh0nk.matplotlib4j.Plot
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.fastq.FastqRecord
import java.io.File
import java.lang.Math
import java.nio.file.Files
import java.nio.file.Paths

fun main(ars: Array<String>) {

    val path_test = "/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first10000_1.fastq.gz"
    val path_run = "/Johnny/data/input/Bacteria/E.coli/K12/ucsd_lane_1/ecoli_mda_lane1.fastq"

    val file = File(path_test)
    val output_dir = "./hw1_output"

    FastqReader(file).use {
        var reads = 0
        val stats = listOf(GcStats(), QualitiesStats(), KmerStats())

        it.iterator().forEach { fastqRecord: FastqRecord ->
            reads++
            val infoCount = 100000
            if (reads % infoCount == 0) {
                println("processed $infoCount reads")
            }

            val readString = fastqRecord.readString
            val quality = fastqRecord.baseQualities

            require(readString.length == quality.size)

            stats.forEach { it.addRead(readString, quality) }
        }

        stats.forEach { it.saveFigures(output_dir) }
    }
}

fun Int.fromPhredScore(): Double = Math.pow(10.0, -this.toDouble() / 10.0)

fun Double.toPhredScore() = Math.round(-10.0 * Math.log10(this))

val qualityCutOff = 4

interface Stats {
    fun addRead(readString: String, quality: ByteArray)
    fun saveFigures(dir: String)
}

class GcStats: Stats {
    private val gcBuckets = 10
    private val gcContents = Array<Int>(gcBuckets + 1, {0})

    override fun addRead(readString: String, quality: ByteArray) {
        var gc = 0
        var nucleotides = 0
        for (i in readString.indices) {
            val nucleotide = readString[i].toLowerCase()
            if (quality[i] < qualityCutOff || readString[i] == 'n') {
                continue
            }
            if (nucleotide == 'g' || nucleotide == 'c') {
                gc++
            }
            nucleotides++

            if (nucleotides != 0) {
                val gcContent = gc.toDouble() / nucleotides
                val gcBucket = Math.floor(gcContent * gcBuckets).toInt()
                gcContents[gcBucket]++
            }
        }
    }

    override fun saveFigures(dir: String) {
        val plt: Plot = Plot.create()

        plt.plot().add(gcContents.indices.map { i -> (i + 0.5) / gcBuckets}, gcContents.map(Int::toDouble))

        Files.createDirectories(Paths.get(dir))
        plt.savefig("$dir/gc.png").dpi(200.0)
        plt.executeSilently()
    }
}

class QualitiesStats: Stats {
    private val readLen = 300

    private val qualitySum: Array<Double> = Array(readLen, { 0.0 })
    private val readCount: Array<Int> = Array(readLen, {0})

    override fun addRead(readString: String, quality: ByteArray) {
        for (i in readString.indices) {
            readCount[i]++
            qualitySum[i] += quality[i].toInt().fromPhredScore()
        }
    }

    override fun saveFigures(dir: String) {
        val plt: Plot = Plot.create()

        val poses: MutableList<Double> = ArrayList()
        val qualities: MutableList<Double> = ArrayList()
        for (i in qualitySum.indices) {
            if (readCount[i] == 0) {
                continue
            }
            poses.add(i.toDouble())
            val phredScore = (qualitySum[i] / readCount[i]).toPhredScore()
            qualities.add(phredScore.toDouble())
        }

        plt.plot().add(poses, qualities)
        plt.savefig("$dir/qual.png")
        plt.executeSilently()
    }
}

class KmerStats: Stats {
    val kmerCount = mutableMapOf<String,Int>()
    override fun addRead(readString: String, quality: ByteArray) {
        for (i in readString.indices) {
            for (k in 1..8) {
                val j = i + k - 1
                if (j >= readString.length || quality[j] < qualityCutOff || readString[j] == 'n') {
                    break
                }
                if (k >= 2) {
                    val kmer = readString.substring(i, i+k)
                    kmerCount[kmer] = (kmerCount[kmer] ?: 0) + 1
                }
            }
        }
    }

    override fun saveFigures(dir: String) {
        for (k in 2..8) {
            val kmerCounts = kmerCount.toList()
                    .filter { (kmer , _) -> kmer.length == k }
                    .sortedByDescending { (_, count) -> count }


            val plt: Plot = Plot.create()
            plt.plot().add(kmerCounts.indices.map(Int::toDouble), kmerCounts.map { (_, count) -> count })
            plt.savefig("$dir/${k}mers.png")
            plt.executeSilently()
        }

    }
}