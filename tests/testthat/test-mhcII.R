test_that("test mhcII prediction works", {
  hu <- mhcII(ag_present = c("DRB1*08:01"),
              ag_stim = c("DQA1_01_01","DQA1_04_01"),
              ag_self = c("DQA1_02_01"),
              seq_len = 15,
              method = "net")

  utx <- "MKSCGVSLATAAAAAAAAAFGDEEKKMAAGKASGESEEASPSLTAEEREALGGLDSRLFGFVRFHEDGARMKALLGKAVRCYESLILKAEGKVESDFFCQLGHFNLLLEDYPKALSAYQRYYSLQSDYWKNAAFLYGLGLVYFHYNAFQWAIKAFQEVLYVDPSFCRAKEIHLRLGLMFKVNTDYESSLKHFQLALVDCNPCTLSNAEIQFHIAHLYETQRKYHSAKEAYEQLLQTENLSAQVKATILQQLGWMHHTVDLLGDKATKESYAIQYLQKSLEADPNSGQSWYFLGRCYSSIGKVQDAFISYRQSIDKSEASADTWCSIGVLYQQQNQPMDALQAYICAVQLDHGHAAAWMDLGTLYESCNQPQDAIKCYLNATRSKNCSNTSGLAARIKYLQAQLCNLPQGSLQNKTKLLPSIEEAWSLPIPAELTSRQGAMNTAQQNTSDNWSGGNAPPPVEQQTHSWCLTPQKLQHLEQLRANRNNLNPAQKLMLEQLESQFVLMQQHQMRQTGVAQVRPTGILNGPTVDSSLPTNSVSGQQPQLPLTRMPSVSQPGVHTACPRQTLANGPFSAGHVPCSTSRTLGSTDTVLIGNNHVTGSGSNGNVPYLQRNAPTLPHNRTNLTSSTEEPWKNQLSNSTQGLHKGPSSHLAGPNGERPLSSTGPSQHLQAAGSGIQNQNGHPTLPSNSVTQGAALNHLSSHTATSGGQQGITLTKESKPSGNTLTVPETSRQTGETPNSTASVEGLPNHVHQVMADAVCSPSHGDSKSPGLLSSDNPQLSALLMGKANNNVGPGTCDKVNNIHPTVHTKTDNSVASSPSSAISTATPSPKSTEQTTTNSVTSLNSPHSGLHTINGEGMEESQSPIKTDLLLVSHRPSPQIIPSMSVSIYPSSAEVLKACRNLGKNGLSNSSILLDKCPPPRPPSSPYPPLPKDKLNPPTPSIYLENKRDAFFPPLHQFCTNPNNPVTVIRGLAGALKLDLGLFSTKTLVEANNEHMVEVRTQLLQPADENWDPTGTKKIWHCESNRSHTTIAKYAQYQASSFQESLREENEKRSHHKDHSDSESTSSDNSGKRRKGPFKTIKFGTNIDLSDDKKWKLQLHELTKLPAFVRVVSAGNLLSHVGHTILGMNTVQLYMKVPGSRTPGHQENNNFCSVNINIGPGDCEWFVVPEGYWGVLNDFCEKNNLNFLMGSWWPNLEDLYEANVPVYRFIQRPGDLVWINAGTVHWVQAIGWCNNIAWNVGPLTACQYKLAVERYEWNKLQNVKSIVPMVHLSWNMARNIKVSDPKLFEMIKYCLLRTLKQCQTLREALIAAGKEIIWHGRTKEEPAHYCSICEVEVFDLLFVTNESNSRKTYIVHCQDCARKTSGNLENFVVLEQYKMEDLMQVYDQFTLAPPLPSASS"
  uty <- "MKSYGLSLTTAALGNEEKKMAAEKARGEGEEGSFSLTVEEKKALCGLDSSFFGFLTRCKDGAKMKTLLNKAIHFYESLIVKAEGKVESDFFCQLGHFNLLLEDYSKALSSYQRYYSLQTDYWKNAAFLYGLGLVYFYYNAFQWAIRAFQEVLYVDPNFCRAKEIHLRLGFMFKMNTDYESSLKHFQLALIDCNVCTLSSVEIQFHIAHLYETQRKYHSAKAAYEQLLQIESLPSQVKATVLQQLGWMHHNMDLIGDNTTKERYAIQYLQKSLEEDPNSGQSWYFLGRCYSCIGKVQDAFVSYRQSIDKSEASADTWCSIGVLYQQQNQPMDALQAYICAVQLDHGHAAAWMDLGILYESCNQPQDAIKCYLNAARSKSCNNTSALTSRIKFLQAQLCNLPQSSLQNKTKLLPSIEEAWSLPIPAELTSRQGAMNTAQQSVSDTWNSVQTASHHSVQQKVYTQCFTAQKLQSFGKDQQPPFQTGSTRYLQAASTNDQNQNGNHTLPQNSKGDAQNHFLRIPTSEEQKIINFTKESKDSRSKSLTSKTSRKDRDTSNICVNAKKHSNHIYQISSVPISSLNNKESVSPDLIIVDNPQLSVLVGETIDNVDHDIGTCDKVNNVHLAIHKKPDNLSASSPSSAISTETLSLKLTEQTHIVTSFISPHSGLHTINGEGHENLESSASVNVGLRPRSQIIPSMSVSIYSSSTEVLKACRSLGKNGLSNGHILLDICPPPRPPTSPYPPLPKEKLNPPTPSIYLENKRDAFFPPLHQFCINPKNPVTVIRGLAGALKLDLGLFSTKTLVEANNEHIVEVRTQLLQPADENWDPSGTKKIWRYENKSSHTTIAKYAQYQACSFQESLREENERRTQVKDYSDNESTCSDNSGRRQKAPFKTIKCGINIDLSDNKKWKLQLHELTKLPAFVRVVSAGNLLSHVGYTILGMNSVQLCMKVPGSRIPGHQENNNFCSVNINIGPGDCEWFVVPEDYWGVLNDFCEKNNLNFLMSSWWPNLEDLYEANVPVYRFIQRPGDLVWINAGTVHWVQAIGWCNNITWNVGPLTAFQYKLAVERYEWNKLQSVKSVVPMVHLSWNMARNIKVSDPKLFEMIKYCLLKILKHCQTLREALVAAGKEVLWHGRINDEPAPYCSICEVEVFNLLFVTNESNSQKTYIVHCQNCARKTSGNLENFVVLEQYKMEDLIQVYDQFTLAPSLSSAS"

  ms <- mhcII(ag_present = c("H2-IAb"),
              ag_self = utx,
              ag_stim = uty,
              seq_len = c("14"),
              nm_self = "utx",
              nm_stim = "uty")
})
