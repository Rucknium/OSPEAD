## Optimal Static Parametric Estimation of Arbitrary Distributions (OSPEAD) for the Monero decoy selection algorithm

This repo contains all submitted documents and code for OSPEAD, most of which has been private until now. OSPEAD is a proposed improvement to Monero's decoy selection algorithm that would reduce the probability that an anti-privacy adversary could correctly guess the real spend in a ring signature.

![Statistical Monero Logo](Statistical-Monero-Logo.gif)

### Results summary

At current Monero ring size of 16, the theoretical minimum attack success through completely random guessing would be 1/16 = 6.25%. According to preliminary estimates, an adversary could take advantage of the divergence between the real spend age distribution and the status quo decoy distribution to achieve an attack success probability of 23.5%, on average, since the August 2022 hard fork. This corresponds to an effective ring size of 4.2. The attack success probability prior to August 2022 may be higher, but this was not measured due to time constraints.

The OSPEAD techniques suggest a new decoy distribution, which would reduce the average attack success probability to 7.6 percent, corresponding to an effective ring size of 13.2.

### Implementation and deployment

It is likely that deployment of a new decoy selection algorithm without a blockchain hard fork would do more harm than good due to some users being slow to upgrade. (For more information about the risk, read my ["Formula for Accuracy of Guessing Monero Real Spends Using Fungibility Defects"](https://github.com/Rucknium/misc-research/blob/main/Monero-Fungibility-Defect-Classifier/pdf/classify-real-spend-with-fungibility-defects.pdf). Therefore, the OSPEAD-derived decoy selection algorithm likely won't be implemented in Monero's standard wallet code before the next hard fork.

Monero's next hard fork is expected to deploy [Full Chain Membership Proofs](https://www.getmonero.org/2024/04/27/fcmps.html), which will eliminate the on-chain ring signature privacy model. However, in certain situations, decoy-based privacy will still be used to provide protection to users' wallets from a potentially malicious spying remote node. Therefore, the OSPEAD-derived decoy distribution can be used in those circumstances. For more details, read ["Initial Probability Density Function for OSPEAD"](CCS-milestone-2/pdf/OSPEAD-milestone-II.pdf).

The OSPEAD documents and code are being publicly released now because there is now an implementable solution to the problems I raised in my original HackerOne submission. Public release will allow greater review and scrutiny of the proposed OSPEAD techniques.

### Documents list

There are three groups of documents, placed in their own directories in the repo:

- `HackerOne-submission` contains my submission to Monero's HackerOne vulnerability response process in September 2021, titled "Research Roadmap for an Overhaul of Monero’s Mixin Selection Algorithm". This was the initial research that suggested a major privacy problem in Monero's decoy selection algorithm, but it did not solve all challenges involving estimating the real spend age distribution, especially the handling of multiple nonstandard decoy selection algorithms being used in the wild. These issues required further research, which was funded by [Monero's Community Crowdfunding System (CCS)](https://ccs.getmonero.org/proposals/Rucknium-OSPEAD-Fortifying-Monero-Against-Statistical-Attack.html).

- `CCS-milestone-1` contains the "Fully Specified Estimation Plan for OSPEAD" document and code submitted in September 2022 as Milestone 1 of the [CCS-funded OSPEAD research project](https://ccs.getmonero.org/proposals/Rucknium-OSPEAD-Fortifying-Monero-Against-Statistical-Attack.html).

- `CCS-milestone-2`contains documents and code submitted in January 2025 as Milestone 2 of the OSPEAD research project. The documents are split into a "Initial Probability Density Function for OSPEAD" PDF and a browser-based website combining procedural narrative and code, viewable at [OSPEAD-docs](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book).

- `old-repo` contains code from the earlier version of this repository.

### How to read the documents

The documents total well over 100 pages. I will suggest portions to read:

- Laypeople (i.e. non-technical people) are suggested to read Sections 19 through 20 (pages 39 to 42) of ["Fully Specified Estimation Plan for OSPEAD"](CCS-milestone-1/pdf/PRIVATE-OSPEAD-Fully-Specified-Estimation-Plan.pdf) and ["Research Roadmap for an Overhaul of Monero’s Mixin Selection Algorithm"](HackerOne-submission/pdf/Roadmap-for-improved-Monero-mixin-selection-algorithm.pdf). Note that the term "mixin" in the latter document changed to "decoy" in the former document.

- Programmers and/or people who prefer to read procedures are suggested to read ["Research Roadmap for an Overhaul of Monero’s Mixin Selection Algorithm"](HackerOne-submission/pdf/Roadmap-for-improved-Monero-mixin-selection-algorithm.pdf) and [OSPEAD-docs](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book). Note that the term for the adversary's attack changed from "Rucknium Ratio Attack" to "Maximum A Postieri (MAP) Decoder attack".

- People interested in the statistical theory underpinning OSPEAD are encouraged to read all of ["Fully Specified Estimation Plan for OSPEAD"](CCS-milestone-1/pdf/PRIVATE-OSPEAD-Fully-Specified-Estimation-Plan.pdf).


### Installing the `decoyanalysis` R package

The `decoyanalysis` R package is a necessary (but not sufficient) condition to reproduce the OSPEAD analysis in Milestone 2. To install `decoyanalysis`, start an R session and input:

```R
install.packages("remotes")
remotes::install_github("Rucknium/OSPEAD",
  subdir = "CCS-milestone-2/decoyanalysis", upgrade = FALSE)
```

On Linux, you may have to install the `liblapack-dev` system package first:

```bash
sudo apt install liblapack-dev
```

### Funding

Most of this research was funded by [Monero's Community Crowdfunding System (CCS)](https://ccs.getmonero.org/proposals/Rucknium-OSPEAD-Fortifying-Monero-Against-Statistical-Attack.html). Thank you to all donors!

## Frequently Asked Questions (FAQs)

**Q: What is the core of the technique that estimates Monero's real spend age distribution?**

**A:** The core of the technique is simple. The real spend distribution is one probability distribution. The decoy distribution is another probability distribution. The on-chain rings are a combination of the decoy and real spend distributions. Given current ring size of 16, the decoy distribution makes up 15 parts of the on-chain data and the real spend makes up 1 part of the on-chain data.  In other words, what is observable on the blockchain is a two-component [mixture distribution](https://en.wikipedia.org/wiki/Mixture_distribution). Using the definition of a mixture distribution, we can simply subtract the decoy distribution from the on-chain data distribution to obtain the real spend distribution. There are other steps more difficult to explain, such as dealing with multiple nonstandard decoy selection algorithms being used in the wild.

**Q: Does OSPEAD base its estimates of Monero's real spend distribution on the data from other cryptocurrencies like Litecoin?**

**A:** No. The OSPEAD estimates only use data from Monero's blockchain and timing information from its transaction pool, which is called a mempool in other blockchains.

In addition to the work directly on the Monero data, I analyzed other blockchains independently to show that the OSPEAD techniques work in realistic scenarios. In a small simulation, I took the Litecoin blockchain data, simulated ring signatures onto it, applied the OSPEAD techniques, and showed that [OSPEAD produced accurate estimates of the original Litecoin distribution](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/successful-simulation.html#plots). I also analyzed other chains to give context about the reasonableness of the OSPEAD estimates of Monero's real spend age distribution, especially the plausibility of large changes in the real spend age distribution from week to week, which are also observed on the transparent UTXO-based payment coins like Bitcoin, Litecoin, Bitcoin Cash, and Dogecoin.

**Q: Does OSPEAD use machine learning techniques?**

**A:** No. OSPEAD uses techniques that are within the traditional statistics paradigm. There is no "training set" or similar idea.

**Q: Does the MAP Decoder attack work by eliminating decoys?**

**A:** No. Decoys aren't completely eliminated. The MAP Decoder attack guesses that the most likely ring member is the real spend. The guess is correct 23.5% of the time.

Think of it this way: there are 16 horses scheduled for a race. The horses are not equally fast. According to the betting markets, one of the horses has a 1-in-4.2 odds of winning. The MAP Decoder attack does _not_ remove 12 out of the 16 horses from the race, and then randomly pick among the remaining 4 with equal probability. Instead, it always bets on the one horse that is most likely to win. It wins the bet (guesses correctly) in 1 out of 4.2 races.

**Q: Are all transactions equally vulnerable to the attack?**

**A:** No. Generally, transactions that spend young outputs are more vulnerable to the attack. Usually, a wallet would spend a young output if a user sends and/or receives XMR more than once in a short period of time. The blue line in Figure 13.5 [here](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/performance-evaluation.html#fig-attack-success-top-1-3) shows the attack success probability according to the age of the real spend, with Monero's current decoy selection algorithm.

The attack success probability also depends on the aggregate behavior of users at any particular time. `OSPEAD-docs` contains [animation of the change in users' spending patterns with respect to output age, for each week](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/real-spend-visualization.html#animations). In some weeks, the average attack success is higher because there is a greater divergence between the decoy selection distribution and the real spend distribution. Sometimes it is lower.

**Q: As a Monero user, how does this affect my privacy?**

**A:** There is no general answer to this. For some people, probabilistic guessing is not a problem. According to certain logic, [plausible deniability](https://en.wikipedia.org/wiki/Plausible_deniability) would still exist even with a high probability of correctly guessing. For other people with extreme threat models, probabalistic methods _are_ a problem for their privacy. Moreover, the MAP Decoder attack could be combined with other attacks to enhance its effectiveness.

**Q: What are the potential known weaknesses of OSPEAD?**

1) One of OSPEAD's steps uses the Bonhomme-Jochmans-Robin (BJR) estimator. The BJR estimator requires that ring members be statistically independent for 100% accuracy of the estimator. Ring members are actually slightly dependent, which has a small negative effect  on the accuracy of the final OSPEAD estimates. The [simulation with Litecoin data](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/successful-simulation.html#plots) shows that the inaccuracy is small. Refer to lines 58-82 on pages 2-3 of ["Initial Probability Density Function for OSPEAD"](https://github.com/Rucknium/OSPEAD/blob/main/CCS-milestone-2/pdf/OSPEAD-milestone-II.pdf) for more explanation.

2) The BJR estimator needs a value for the number of distinct decoy selection algorithms being used in the wild. This was difficult to obtain directly from the on-chain data, so a guess of 4 was made. The [simulation with Litecoin data](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/successful-simulation.html#plots) shows that OSPEAD still performs well if the number of distinct algorithms is overestimated. However, underestimation could create inaccuracies.

3) In 2023, [an off-by-one-block bug was discovered and patched in the standard implementation of Monero's decoy selection algorithm](https://www.getmonero.org/2023/06/08/10block-old-decoy-selection-bug.html). It is unknown how many user have updated their wallets to the patched version. Therefore, there are two "standard" decoy selection algorithms operating in the wild. The difference is only one block (about two minutes on average), but it affects the youngest outputs that are the most likely to be spent. I made a reasonable guess of the adoption pace of the patched wallet version. The reasonable guess may be changed in a future, final version of OSPEAD. If the guess is very wrong, then the OSPEAD estimate of the real spend distribution may be inaccurate in the first two minutes of the distribution. See [here](https://rucknium.github.io/OSPEAD/CCS-milestone-2/OSPEAD-docs/_book/transformed-age.html) for more information.

4) The OSPEAD statistical code itself needs to be triple-checked for off-by-one errors.


