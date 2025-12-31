library(data.table)
library(ggplot2)
#devtools::install_github("frhl/bravastring")
library(bravastring)
options(repr.matrix.max.cols=50, repr.matrix.max.rows=10)

l <- lapply(files[2:16], function(f){
        d <- fread(f)
        d[ , heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
        d[ , id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
        d[, annotation := stringr::str_extract(f, "(?<=continuous_\\d\\.).+?(?=\\.combined)")]
        # restrict to variants with at least one biallelic variant
        d <- d[!is.na(d$REC.P),]
        return(d)
    })

d[, REC.PEXPT := get_expected_p(REC.P), by = .(heritability, id)]


ggplot(d, aes(x=-log10(REC.PEXPT), y=-log10(REC.P))) + 
    geom_point() +
    geom_abline(linetype='dashed') +
    facet_wrap(~heritability+id)


# get colors
colors <- list(red=c("#E35278"), orange=c("#E2A98C"), green=c("#9EC0A6"), blue=c("#009894"))

bic_color_scale <- scale_color_manual(
    values = c(
    "<5" = colors$red,
    "[5, 10)" = colors$orange,
    "[10, 25)" = colors$green,
    ">25" = colors$blue
  ))
#labels = c("<0.05", "[0.05, 1)","[1, 10)", ">10")
ebc_color_scale <- scale_color_manual(
    values = c(
    "<0.02" = colors$red,
    "[0.02, 1)" = colors$orange,
    "[1, 10)" = colors$green,
    ">10" = colors$blue
  ))
  
  set.seed(152)
idx <- sample(1:nrow(d), 10000, replace=FALSE)
d_plot <- d[idx,]
d_plot <- d
#d_plot[, p.value.expt := get_expected_p(p.value), by = .(BIC_interval, prevalence, test_type)]
d_plot[, p.value.expt := get_expected_p(p.value), by = .(EBC_interval, prevalence, test_type)]

# define ribbon that can exits seperately
n <- max(aggregate(p.value~prevalence+test_type+BIC_interval, FUN=length, data=d_plot)$p.value)
seq.p.value <- seq(0, 1, length.out=n)
seq.p.value.expt <- get_expected_p(seq.p.value)
dt_ribbon <- data.table(
    p.value = -log10(seq.p.value),
    p.value.expected = -log10(seq.p.value.expt),
    clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n:1, shape1 = 1:n)),
    cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n:1, shape1 = 1:n))
)

options(repr.plot.width=12, repr.plot.height=8)
p <- ggplot(d_plot,
        aes(x=-log10(p.value.expt),
            y=-log10(p.value), col=EBC_interval)
        ) +
    geom_ribbon(
        data=dt_ribbon,
        aes(x=p.value.expected, y=p.valuOBOBe, ymax=cupper, ymin=clower),
        fill="grey90", color="black") +
    ebc_color_scale +
    geom_point_rast(size=3) +
    facet_wrap(~prevalence) +
    geom_abline(intercept=0, slope=1, color='black', linetype = "dashed") +
    xlab(TeX("$-\\log_{10}(P_{observed})$")) +
    ylab(TeX("$-\\log_{10}(P_{expected})$")) +
    theme_classic() +
    theme(
        legend.position = "right",
        strip.text = element_text(size=16),
        axis.text=element_text(size=13),
        axis.title=element_text(size=16,face="bold"),
        axis.title.x = element_text(margin=ggplot2::margin(t=15)),
        axis.title.y = element_text(margin=ggplot2::margin(r=15)),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
    )
