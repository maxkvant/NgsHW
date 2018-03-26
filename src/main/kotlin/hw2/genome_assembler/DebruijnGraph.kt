package hw2.genome_assembler

import K
import Kmer
import reverseComplement

class DebruijnGraph(kmers1: Set<Kmer>) {
    private val kmers: Array<Kmer> = kmers1.toTypedArray()
    private val kmerId: MutableMap<Kmer, Int> = hashMapOf()
    private val nucleotidesSts = arrayOf("a", "c", "g", "t")
    private val vertices: MutableSet<Int> = kmers.indices.toMutableSet()
    private val revVertex: IntArray = IntArray(kmers.size)

    private val edges: Array<Array<Edge?>> = Array(kmers.size, { Array<Edge?>(nucleotidesSts.size, { null }) })

    init {
        for (i in kmers.indices) {
            require(kmers[i].length == K)
            kmerId[kmers[i]] = i
        }
        for (i in kmers.indices) {
            revVertex[i] = kmerId[reverseComplement(kmers[i])]!!
        }
    }


    fun addEdge(edgeStr: String) {
        require(edgeStr.length == K + 1)
        if (edgeStr.contains('n')) {
            return
        }

        val c = nucleotideCode(edgeStr[K])

        val from = kmerId[edgeStr.substring(0, K)]!!
        val to = kmerId[edgeStr.substring(1, K + 1)]!!

        if (edges[from][c] == null) {
            edges[from][c] = Edge(from, to, nucleotidesSts[c], 1)
        }

        edges[from][c]!!.coverage += 1
    }

    fun removeEdges(shouldRemove: (Edge) -> Boolean) {
        for (v in kmers.indices) {
            for (c in nucleotidesSts.indices) {
                val edge = edges[v][c]
                if (edge != null && shouldRemove(edge)) {
                    removeEdge(edge)
                }
            }
        }
        removeUnusedVertices()
    }

    fun removeTails() {
        val curEdges = allEdges()
        val longestCoverage = curEdges.maxBy { it.str.length }!!.coverageAverage
        val lenCut = 1.5 * K

        removeEdges { edge ->
            (degIn(edge.from) == 0 || degOut(edge.to) == 0)
            && (edge.coverageAverage < longestCoverage / 2.0)
            && (edge.str.length < lenCut)
        }

        contract()

        removeEdges { edge ->
            (degIn(edge.from) == 0 || degOut(edge.to) == 0)
            && (edge.str.length < lenCut)
        }

        contract()
    }

    private fun removeEdge(edge: Edge) {
        val revE = edge.revEdge()
        edges[revE.from][revE.firstNucleotideCode] = null
        edges[edge.from][edge.firstNucleotideCode] = null
    }

    fun contract() {
        for (v in kmers.indices) {
            if (!vertices.contains(v) || degOut(v) < 2) {
                continue
            }

            contractPathsFrom(v)
        }

        for (v in kmers.indices) {
            if (!vertices.contains(v) || degOut(v) != 1) {
                continue
            }

            contractPathsFrom(v)
        }

        removeUnusedVertices()
    }

    fun allEdges(): List<Edge> = vertices.flatMap { v -> edges[v].filterNotNull() }

    private fun degOut(v: Int): Int {
        return edges[v].count { it != null}
    }

    private fun degIn(v: Int): Int {
        return degOut(revVertex[v])
    }

    private fun nucleotideCode(nucleotide: Char): Int =
            when (nucleotide) {
                'a' -> 0
                'c' -> 1
                'g' -> 2
                't' -> 3
                else -> throw IllegalArgumentException("$nucleotide isn't nucleotide")
            }

    private fun removeUnusedVertices() {
        for (v in kmers.indices) {
            if (degIn(v) == 0 && degOut(v) == 0) {
                vertices.remove(v)
            }
            if (!vertices.contains(v)) {
                edges[v].fill(null)
            }
        }
    }

    private fun contractPathsFrom(v: Int) {
        for (c in edges[v].indices) {
            val edge = edges[v][c]
            if (edge != null) {
                val used: MutableSet<Int> = mutableSetOf(v)

                val sb = StringBuilder()
                sb.append(edge.str)

                var u = edge.to
                var coverage: Long = edge.coverage

                while (degIn(u) == 1 && degOut(u) == 1 && !used.contains(u)) {
                    used.add(u)
                    vertices.remove(u)

                    val curEdge: Edge = edges[u].first { it != null }!!
                    coverage += curEdge.coverage
                    sb.append(curEdge.str)

                    u = curEdge.to
                }

                val str = sb.toString()
                val newEdge = Edge(v, u, str, coverage)
                edges[v][c] = newEdge

                vertices.add(v)
                vertices.add(u)
            }
        }
    }

    inner class Edge(val from: Int, val to: Int, val str: String, var coverage: Long) {
        val coverageAverage: Double
            get() = coverage.toDouble() / str.length

        val fromStr: Kmer
            get() = kmers[from]

        val toStr: Kmer
            get() = kmers[to]

        val fullString: String
            get() = kmers[from] + str

        fun revEdge(): Edge {
            val newFrom = revVertex[to]
            val newTo = revVertex[from]
            for (edge in edges[newFrom]) {
                if (edge != null && edge.to == newTo) {
                    return edge
                }
            }
            throw IllegalStateException("rev edge not found")
        }

        internal val firstNucleotideCode: Int
            get() = nucleotideCode(str[0])
    }
}