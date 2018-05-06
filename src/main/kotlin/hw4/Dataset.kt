package hw4

data class Dataset(val reference: String, val reads: List<String>) {
    init {
        require(reads.size == 1 || reads.size == 2)
    }
}