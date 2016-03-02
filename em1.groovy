#!/usr/bin/env groovy

/*
 * A 1-locus EM implementation to estimate haplotype frequencies.
 *
 * usage: em1.groovy <input file> <frequencies output> <individual output>
 *   e.g., em1.groovy sout_interps.txt output/em-haps.txt output/em-ids.txt
 *
 * The input file is a tab-separted 3-column file. Each line
 * represents an individual's haplotype list.
 *   column 1: an individual ID
 *   column 2: not currently used; reserved for allelic information
 *   column 3: a haplotype list (e.g., 1+1|2+3); all ambiguity is in the '|'
 *
 * The frequencies output is a tab-separated 2-column file. Each line
 * contains the haplotype name and its estimated frequency.
 *
 * The individual output is a tab-separated 3-column file. Each line
 * contains the individual ID, its estimated haplotype-pair list (which may
 * be ambiguous), and a single haplotype pair, which is the same as
 * column 2 if the haplotype list is unambiguous or a randomly chosen
 * haplotype pair from column 2 if the haplotype list is ambiguous.
 *
 * @author Dave Roe
 * @author Rui Kuang
 */

import java.util.Random  

// things that may change per run
debugging = false

// misc vars
err = System.err
if((args.size() == 0) || (args[0] == "-h")) {
    err.println "usage: ./em1.groovy <input file> <frequencies output> <individual output>"
    System.exit(1)
}       
inputFile = args[0]
hapOutFile = args[1]
idOutFile = args[2]

// load input
(idList, hlList) = loadHapPairLists(inputFile)
if(idList.size() == 0) {
    System.exit(1)
}

// h: a frequency map
// to start, each haplotype gets an equal expected frequency
Map<String, Float> h = [:]
aSum = 1
eDelta = true
round = 0
while(eDelta == true) {
    round++;
    h = m(hlList, h, aSum)                // M-step
    (hlList, eDelta, aSum) = e(hlList, h) // E-step
}  // while eDelta
err.println "done in ${round} rounds"

// output the hap predictions per individual
err.println "outputting ${idOutFile}..."
idOut = new PrintWriter(new File(idOutFile).newOutputStream(), true)
idOut.println "id\thaplotype list\thaplotype pair"
hlList.eachWithIndex { hl, i ->
    id = idList[i]
    randHP = randomHapPair(hl)
    idOut.println "${id}\t${hl}\t${randHP}"
} // each 

// output the haplotype frequency estimations
err.println "done. outputting ${hapOutFile}..."
hapOut = new PrintWriter(new File(hapOutFile).newOutputStream(), true)
hapOut.println "haplotype\tfrequency"
h.each { hap, freq ->
    hapOut.println "${hap}\t${freq.round(4)}"
}
err.println "done"
// end main

/*
 * loadHapPairLists
 * Loads the input into two lists.
 *
 * @param inputFileName full path to file; The input file is a tab-separted 
 *                      3-column file. Each line represents an individual's 
 *                      haplotype list.
 *                      column 1: an individual ID
 *                      column 2: not used
 *                      column 3: a haplotype list (e.g., 1+1|2+3); 
 *                                all ambiguity is in the '|'
 * @return a List with two Lists; the first contains the IDs, the second contains
 *         the haplotype pair lists in the same order as the IDs
 */
def ArrayList<ArrayList> loadHapPairLists(String inputFileName) {
    err.println "loading file $inputFile ..."
    in_csv = new FileReader(inputFile)
    // e.g., haplotype list: 1+1|2+3
    ArrayList<String> idList = []
    ArrayList<String> hlList = []

    c=0;
    in_csv.eachLine { row ->
        c++;
        (id, alist, slist) = row.split()
        if(debugging) { 
            err.println "${id}"
        }

        idList.add(id)
        hlList.add(slist)
    } // each input line
    err.println "done: ${idList.size()} IDs"

    return [idList, hlList]
} // loadHapPairLists

/*
 * M-step
 *
 * @param hlList List of haplotype-pair lists (e.g., 1+1|2+3)
 * @param h Map of haplotype names to frequencies
 * @param aSum sum of highest scores for each individual (a_ijk)
 * @return h
 */
def Map m(ArrayList hlList, Map h, Double aSum) {
    if(debugging) {
        err.println "entering M-step..."
    }
    Map<String, Integer> sumMap = [:] // haplotype name -> sum
    sumMapTotal = 0.0
    hlList.eachWithIndex { hl, i ->      // each individual
        if(debugging) {
            err.println "${idList[i]}=${hl}"
        }
        hpList = hl.split('\\|')
        // even shares per haplotype pair
        hpListSize =hpList.size()
        share = 1/(hpListSize*2)
        Double a = 1
        if(debugging) {
            err.println "share=${share}, a=${a}"
        }
        hpList.each { hp -> // each haplotype pair
            hapList = hp.split('\\+')
            hapListSize = hapList.size()
            hapList.each {  hap ->
                thisA = a*share
                if(sumMap[hap] == null) {
                    sumMap[hap] = thisA
                } else {
                    sumMap[hap] += thisA
                }
                if(debugging) {
                    err.println "m: thisA=${thisA}, sumMap[${hap}]=${sumMap[hap]}"
                }
                sumMapTotal += thisA
            } // each haplotype in a pair
        } // each haplotype pair
    } // each individual's haplotype list

    // average the a sums per haplotype across all haplotypes
    Float totalSum = 0.000
    h.keySet().each { hap ->
        h[hap] = (Double)0.0
    }
    // normalize
    k = sumMap.keySet()
    k.each { hap ->
        if(debugging) {
            err.println "hap=${hap}, h[${hap}]=${h[hap]}"
        }
        num = sumMap[hap]
        r = num/sumMapTotal
        if(debugging) {
            err.println "new h[${hap}]=" + r + ", ${num}/${sumMapTotal}"
        }
        h[hap] = r
    }
    if(debugging) {
        err.println "sumMapTotal=${sumMapTotal}"
        err.println "done: M-step"
    }
    return h
} // m

/*
 * E-step
 *
 * @param hlList List of haplotype-pair lists (e.g., 1+1|2+3)
 * @param h Map of haplotype names to frequencies
 * @return a List with four items
 *         hlList 
 *         Boolean: true if hlList changed; false if not
 *         a List of highest scores (a_ijk) per individual; each a_ijk
 *           should be used for _both_ haplotypes
 *         aSum: sum of highest scores for each individual (a_ijk)
 */
def ArrayList e(ArrayList hlList, Map h) {
    if(debugging) {
        err.println "entering E-step..."
    }
    Double aSum = 0
    eDelta = false // did any frequency change during the E-step?
    hlListNew = []
    Double maxFreq   // frequency of the most probably pair of haplotypes
    hlList.each { hl ->      // each individual
        if(debugging) {
            err.println "e: hl=${hl}"
        }
        maxFreq = 0
        hpList = hl.split('\\|')
        maxFreqList = [] // list of most likely haplotype pairs
        hpList.each { hp -> // each haplotype pair
            (h1, h2) = hp.split('\\+')
            if(debugging) {
                err.println "e: h1=${h1}, h2=${h2}"
            }
            freq = h[h1] * h[h2]
            if(h1 != h2) {
                freq = freq * 2
            }
            if(freq > maxFreq) { // delete the old max list, if any
                if(debugging) {
                    err.println "e: max now ${h1}+${h2}; ${h1}(${h[h1]}) * ${h2}(${h[h2]}) = ${freq}"
                }
                maxFreqList = []
            }
            if(freq >= maxFreq) {
                maxFreq = freq
                maxFreqList.add(hp)
                if(debugging) {
                    err.println "max: ${h1}(${h[h1]}) * ${h2}(${h[h2]}) = ${freq}"
                }
            } else { // this haplotype pair is not highest scoring
                eDelta = true
            }
        } // each haplotype pair

        aSum += maxFreq // this is for both haplotypes
        listMF = maxFreqList.join('|')
        if(debugging) {
            err.println "m: maxFreq=${maxFreq}, listMF=${listMF}"
        }
        hlListNew.add(listMF)
    } // each individual

    if(debugging) {
        err.println "done: E-step, aSum=${aSum}"
    }

    hlList = hlListNew
    return [hlList, eDelta, aSum]
} // e

/*
 * randomHapPair
 * Returns a random haplotype pair from a haplotype list.
 *
 * @param hl a haplotype list
 * @ret a randomly-chosen single pair of haplotypes from the haplotype list
 */
def String randomHapPair(String hl) {
    hpList = hl.split('\\|')
    size = hpList.size()
    if(size > 1) {
        Random rand = new Random()
        i = rand.nextInt(size)
        hl = hpList[i]
    }
    return hl
} // randomHapPair
