
typealias Kmer = String
val K = 55

fun reverseComplement(str: String): String {
    val chars = str.toCharArray()
    chars.reverse()
    for (i in chars.indices) {
        val c = chars[i]
        chars[i] = when (c) {
            'a' -> 't'
            't' -> 'a'
            'g' -> 'c'
            'c' -> 'g'
            else -> 'n'
        }
    }
    return String(chars)
}