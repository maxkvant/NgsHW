package common

fun execCmd(cmd: String): Int {
    println(cmd)
    val p = Runtime.getRuntime().exec(arrayOf("/bin/sh", "-c", cmd))
    val exitCode: Int = p.waitFor()
    p.getInputStream().bufferedReader().use { reader ->
        reader.lineSequence().forEach { println(it) }
    }
    return exitCode
}