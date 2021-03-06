package hw1

import com.github.sh0nk.matplotlib4j.Plot
import htsjdk.samtools.fastq.FastqReader
import htsjdk.samtools.fastq.FastqRecord
import java.io.File
import java.lang.Math
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.math.roundToInt

interface Stats {
    //добавить инфомацию о риде (для подсчёта статистики)
    fun addRead(readString: String, quality: ByteArray)

    //сохранить графики в папку dir
    fun saveFigures(dir: String)
}

fun main(ars: Array<String>) {
    val pathRun = "/Johnny/data/input/Bacteria/E.coli/K12/ucsd_lane_1/ecoli_mda_lane1.fastq"
    //val pathTest = "/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first10000_1.fastq.gz"

    val file = File(pathRun)
    val outputDir = "./hw1_output"
    val oldIllumina = true

    FastqReader(file).use {
        var reads = 0
        val stats: List<Stats> = listOf(GcStats(), QualityStats(), KmerSpectrumStats())

        it.iterator().forEach { fastqRecord: FastqRecord ->
            reads++
            val infoCount = 100000
            if (reads % infoCount == 0) {
                println("processed $reads reads")
            }

            val readString = fastqRecord.readString
            val quality = fastqRecord.baseQualities
            if (oldIllumina) {
                for (i in quality.indices) {
                    quality[i] = (quality[i] - (64 - 33)).toByte()
                }
            }

            require(readString.length == quality.size)

            stats.forEach { it.addRead(readString, quality) }
        }

        Files.createDirectories(Paths.get(outputDir))
        stats.forEach { it.saveFigures(outputDir) }
    }
}

fun Int.fromPhredScore(): Double = Math.pow(10.0, -this.toDouble() / 10.0)

fun Double.toPhredScore(): Double = -10.0 * Math.log10(this)

val qualityCutOff = 4

//строит GC график
class GcStats: Stats {
    private val gcBuckets = 100
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
        }
        if (nucleotides > readString.length * 0.5) {
            val gcBucket = gc * gcBuckets / nucleotides
            gcContents[gcBucket]++
        }
    }

    override fun saveFigures(dir: String) {
        val plt: Plot = Plot.create()
        plt.ylabel("# reads")
        plt.xlabel("% gc")
        plt.plot().add(gcContents.indices.map { i -> (i + 0.5) / gcBuckets},
                       gcContents.toList(),
                "r+")

        plt.savefig("$dir/gc.png").dpi(200.0)
        plt.executeSilently()
    }
}

//Строит графики среднего качества, средней вероятности ошибки
class QualityStats : Stats {
    private val readLen = 300

    private val qualitySum: Array<Double> = Array(readLen, { 0.0 })
    private val probabilitySum: Array<Double> = Array(readLen, { 0.0 })
    private val readCount: Array<Int> = Array(readLen, {0})

    override fun addRead(readString: String, quality: ByteArray) {
        for (i in readString.indices) {
            readCount[i]++
            qualitySum[i] += quality[i].toDouble()
            probabilitySum[i] += quality[i].toInt().fromPhredScore()
        }
    }

    override fun saveFigures(dir: String) {
        val plt: Plot = Plot.create()

        val poses: MutableList<Double> = ArrayList()
        val qualities: MutableList<Double> = ArrayList()
        val qualities2: MutableList<Double> = ArrayList()
        val probabilities: MutableList<Double> = ArrayList()
        for (i in qualitySum.indices) {
            if (readCount[i] == 0) {
                continue
            }
            poses.add(i.toDouble())
            qualities.add(qualitySum[i] / readCount[i])
            val pAverage: Double = probabilitySum[i] / readCount[i]
            qualities2.add(pAverage.toPhredScore())
            probabilities.add(pAverage)
        }


        plt.xlabel("position")
        plt.ylabel("quality average")
        plt.ylim(0, qualities.max()!! + 1)
        plt.plot().add(poses, qualities, ".").label("quality average")
        plt.legend().loc(0)
        plt.plot().add(poses, qualities2, "r+").label("phred_score(probability_average)")
        plt.legend().loc(  0)
        plt.savefig("$dir/qual_q.png")
        plt.executeSilently()
        plt.close()

        plt.xlabel("position")
        plt.ylabel("probability average")
        plt.ylim(0, probabilities.max()!! * 1.03)
        plt.plot().add(poses, probabilities, "v")
        plt.savefig("$dir/qual_p.png")
        plt.executeSilently()
    }
}

//Строит спектр k-меров
class KmerSpectrumStats : Stats {
    val kmerCount = Array(9, { i ->
        val size = Math.pow(5.0, i.toDouble()).roundToInt() + 1
        Array<Int>(size, {0})
    })

    override fun addRead(readString: String, quality: ByteArray) {
        for (i in readString.indices) {
            var hash = 0
            for (k in 1..8) {
                val j = i + k - 1
                if (j >= readString.length) {
                    break
                }
                val nucleotide = readString[j].toLowerCase()
                if (quality[j] < qualityCutOff || nucleotide == 'n') {
                    break
                }

                val nucleotideCode = when (nucleotide) {
                    'a' -> 1
                    'c' -> 2
                    'g' -> 3
                    't' -> 4
                    else -> 0
                }
                hash = hash * 5 + nucleotideCode

                kmerCount[k][hash]++
            }
        }
    }

    override fun saveFigures(dir: String) {
        for (k in 2..8) {
            val kmerCounts = kmerCount[k].toList()
                .filter { it > 0 }
            val maxCount = kmerCounts.max() ?: 0
            val buckets = 200
            val distribution = Array(buckets + 1, {0.0})
            kmerCounts.forEach {
                val bucket = (it.toLong() * buckets) / maxCount
                distribution[bucket.toInt()] += 1.0 / kmerCounts.size * buckets
            }

            val plt: Plot = Plot.create()
            plt.title("${k}mer spectrum")
            plt.ylabel("density")
            plt.xlabel("# copy number")
            plt.ylim(0, distribution.max()!! + 1)
            plt.plot().add(
                    distribution.indices.map { (it + 0.5) * maxCount / buckets },
                    distribution.toList(),
                    "r+"
            )
            plt.savefig("$dir/${k}mers.png")
            plt.executeSilently()
        }

    }
}