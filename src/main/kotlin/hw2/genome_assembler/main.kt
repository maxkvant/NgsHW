package hw2.genome_assembler

import K
import Kmer
import guru.nidi.graphviz.engine.Format
import guru.nidi.graphviz.engine.Graphviz
import guru.nidi.graphviz.parse.Parser
import java.io.File
import kotlin.math.*

val outDir = "hw2_output"

fun buildGraph(reader: SeqReader): DebruijnGraph {
    val kmers: MutableSet<Kmer> = mutableSetOf()

    reader.forCharArraysComplement { chars ->
        for (i in 0 .. chars.size - K) {
            val kmer: Kmer = String(chars, i, K)
            if (!kmer.contains('n')) {
                kmers.add(kmer)
            }
        }
    }

    val debruijnGraph = DebruijnGraph(kmers)

    reader.forCharArraysComplement { chars ->
        for (i in 0 .. chars.size - (K + 1)) {
            val edgeStr = String(chars, i, K + 1)
            debruijnGraph.addEdge(edgeStr)
        }
    }

    val lowCoverage = min(10.0, 0.1 * debruijnGraph.allEdges().map { it.coverageAverage }.average())
    debruijnGraph.removeEdges {
        it.coverageAverage < lowCoverage
    }
    debruijnGraph.contract()
    debruijnGraph.removeTails()
    return debruijnGraph
}

fun outputPng(dotStr: String, fileName: String) {
    val gr = Parser.read(dotStr)
    val outfile = File("$outDir/$fileName.png")
    val dotFile = File("$outDir/$fileName.dot")
    println(Format.PNG)
    dotFile.writeText(dotStr)
    Graphviz.fromGraph(gr).render(Format.PNG).toFile(outfile)
}

fun outputFasta(edges: List<DebruijnGraph.Edge>, fileName: String) {
    val outfile = File("$outDir/$fileName.fasta")
    val sb = StringBuilder()
    edges.sortedByDescending { it.str.length }
            .withIndex()
            .forEach { (i, edge) ->
                sb.appendln(">edge_${i}_len_${edge.fullString.length}_cov_${edge.coverageAverage} ${edge.fromStr} -> ${edge.toStr}")
                sb.appendln(edge.fullString)
                sb.appendln()
            }
    outfile.writeText(sb.toString())

}

fun toDotStr(edges: List<DebruijnGraph.Edge>): String {
    val vertices = edges.flatMap { listOf(it.fromStr, it.toStr) } .toSet()
    val sb = StringBuilder()
    vertices.forEach {
        sb.appendln(it)
    }
    sb.appendln()
    edges.forEach { edge ->
        val label = "len=${edge.str.length}\ncov=${edge.coverageAverage.roundToInt()}"
        sb.appendln("${edge.fromStr} -> ${edge.toStr} [label=\"$label\"]")
    }
    return "digraph debruijn {\n" + sb.toString() + "}\n"
}

fun main(args: Array<String>) {
    val inputDir = "/media/maxim/DATA/Downloads/NGS/data"

    val pathsAll = listOf(
    //        "ECOLI_IS220_QUAKE_1K_paired_reads.fasta"
    //        , "ECOLI_IS220_QUAKE_1K_single_reads.fasta"
            "s_6.first10000.fastq"
            , "s_6.first1000.fastq"
    //        , "test1.fasta"
    //        , "test2.fasta"
            , "s_6.first100000.fastq"
    )


    pathsAll.forEach { path ->
        val file = File("$inputDir/$path")
        val name = path.substring(0, path.length - ".fast_".length)
        val reader = when {
            path.endsWith(".fasta") -> FastaSeqReader(file)
            path.endsWith(".fastq") -> FastqSeqReader(file)
            else -> throw RuntimeException()
        }

        val debruijnGraph = buildGraph(reader)
        val edges = debruijnGraph.allEdges()

        outputFasta(edges, name)
        outputPng(toDotStr(edges), name)
    }
}