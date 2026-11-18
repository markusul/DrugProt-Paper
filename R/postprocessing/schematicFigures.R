# Figures 1 and 5 are schematic drawings.
# Figure 1 shows three screenshots of the application from schematics/, each with an icon, and a sketch of the model.
# Figure 5 is the causal graph of the Methods.
source("R/figureUtils.R")

##### Figure 1: overview of the software and of the model #####
# all coordinates are in inches, panel (a) is 3.2 in wide and panel (b) 4.3 in

# polygon of a rectangle with rounded corners
round_rect <- function(x0, x1, y0, y1, r, n = 12) {
  a <- seq(0, pi / 2, length.out = n)
  data.frame(x = c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a)),
             y = c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a)))
}

# a monitor that shows a screenshot, centered at (x, y)
monitor <- function(file, x, y, w = 0.95, h = 0.6) {
  img <- png::readPNG(file)
  scale <- min(w / ncol(img), h / nrow(img))
  iw <- ncol(img) * scale
  ih <- nrow(img) * scale
  y <- y + 0.07
  list(
    geom_polygon(data = round_rect(x - w / 2 - 0.035, x + w / 2 + 0.035, y - h / 2 - 0.035, y + h / 2 + 0.035, 0.04),
                 aes(x, y), fill = "grey20"),
    annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2, fill = "white"),
    annotation_raster(img, xmin = x - iw / 2, xmax = x + iw / 2, ymin = y - ih / 2, ymax = y + ih / 2,
                      interpolate = TRUE),
    annotate("rect", xmin = x - 0.06, xmax = x + 0.06, ymin = y - h / 2 - 0.13, ymax = y - h / 2 - 0.03, fill = "grey20"),
    geom_polygon(data = round_rect(x - 0.2, x + 0.2, y - h / 2 - 0.17, y - h / 2 - 0.125, 0.02), aes(x, y), fill = "grey20")
  )
}

fig_arrow <- function(x, y, xend, yend, linewidth = 0.6, length = 0.08) {
  annotate("segment", x = x, y = y, xend = xend, yend = yend, linewidth = linewidth,
           arrow = arrow(length = unit(length, "in"), type = "closed"))
}
fig_text <- function(x, y, label, size = 9, ...) {
  annotate("text", x = x, y = y, label = label, family = fig_family, size = size / .pt, lineheight = 0.95, ...)
}

# a protein as a small ribbon diagram, a sheet of three strands between two helices,
# of width about 2 * s, centered at (x, y) and rotated by angle degrees
# halo = TRUE puts a white disc behind the protein, which hides the lines that lead to it
protein <- function(x, y, s = 0.26, linewidth = 0.45, angle = -30, halo = FALSE, colour = "black") {
  # a Bezier curve from p0 to p3 with the control points p1 and p2
  bezier <- function(p0, p1, p2, p3) {
    u <- seq(0, 1, length.out = 40)
    b <- cbind((1 - u)^3, 3 * (1 - u)^2 * u, 3 * (1 - u) * u^2, u^3) %*% rbind(p0, p1, p2, p3)
    data.frame(x = b[, 1], y = b[, 2])
  }
  # one half of the chain: a tail, a helix, and a loop to the lower end of a strand at height end
  t <- seq(pi, 9 * pi, length.out = 200)
  half <- function(end) {
    rbind(bezier(c(-0.80, 0.95), c(-0.94, 0.90), c(-0.99, 0.72), c(-0.99, 0.5)),
          data.frame(x = -0.82 + 0.17 * cos(t), y = 0.5 - (t - pi) / (8 * pi) + 0.11 * sin(t)),
          bezier(c(-0.99, -0.5), c(-0.99, -0.84), c(-0.36, end - 0.34), c(-0.36, end)))
  }
  # the second half is the first half turned by 180 degrees
  chain <- rbind(cbind(half(-0.55), group = 1), cbind(-half(-0.62), group = 2),
                 # the two turns between the strands
                 cbind(bezier(c(-0.36, 0.62), c(-0.36, 0.88), c(0, 0.88), c(0, 0.55)), group = 3),
                 cbind(bezier(c(0, -0.62), c(0, -0.88), c(0.36, -0.88), c(0.36, -0.55)), group = 4))
  # the three strands as arrows that point up, down, and up
  arrow_shape <- data.frame(x = c(-0.08, -0.08, -0.16, 0, 0.16, 0.08, 0.08),
                            y = c(-0.55, 0.34, 0.34, 0.62, 0.34, 0.34, -0.55))
  strands <- do.call(rbind, lapply(1:3, function(k) {
    data.frame(x = arrow_shape$x + 0.36 * (k - 2), y = arrow_shape$y * (-1)^(k + 1), group = k)
  }))
  place <- function(d) {
    a <- angle * pi / 180
    data.frame(x = x + s * (cos(a) * d$x - sin(a) * d$y), y = y + s * (sin(a) * d$x + cos(a) * d$y), group = d$group)
  }
  a <- seq(0, 2 * pi, length.out = 60)
  disc <- data.frame(x = x + 1.25 * s * cos(a), y = y + 1.25 * s * sin(a))
  list(if (halo) geom_polygon(data = disc, aes(x, y), fill = "white"),
       geom_path(data = place(chain), aes(x, y, group = group), colour = colour, linewidth = linewidth,
                 lineend = "round", linejoin = "round"),
       geom_polygon(data = place(strands), aes(x, y, group = group), fill = "white", colour = colour,
                    linewidth = linewidth, linejoin = "mitre"))
}

# a drug capsule, half filled, rotated by the given angle
capsule <- function(x, y, angle = pi / 4, fill = "black", len = 0.22, r = 0.1, linewidth = 0.6, n = 20) {
  a <- seq(pi / 2, 3 * pi / 2, length.out = n)
  left <- cbind(-len / 2 + r * cos(a), r * sin(a))
  outline <- rbind(left, cbind(len / 2 - r * cos(a), -r * sin(a)))
  half <- rbind(left, c(0, -r), c(0, r))
  rot <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
  shape <- function(p) data.frame(x = x + (p %*% rot)[, 1], y = y + (p %*% rot)[, 2])
  list(geom_polygon(data = shape(outline), aes(x, y), fill = "white", colour = "black", linewidth = linewidth),
       geom_polygon(data = shape(half), aes(x, y), fill = fill, colour = "black", linewidth = linewidth))
}

# (a) the three steps of a query, each with a screenshot and an icon
rows <- c(2.55, 1.55, 0.55)
net <- expand.grid(k = 1:3, t = 1:3)
net$x <- 3.07 + c(-0.26, 0, 0.26)[net$t]
net$y <- rows[3] + c(0.26, 0, -0.26)[net$k]
net_links <- data.frame(from = c(1, 1, 2, 3, 3, 4, 5, 5, 6), to = c(4, 5, 6, 5, 6, 8, 7, 9, 9))

fig1_a <- ggplot() +
  fig_text(0.72, rows[1], "Select proteins\nof interest") +
  fig_text(0.72, rows[2], "1. Drug effects\nand interactions") +
  fig_text(0.72, rows[3], "2. Temporal\nprotein network") +
  monitor("schematics/screen_select.png", 2.0, rows[1]) +
  monitor("schematics/screen_drugs.png", 2.0, rows[2]) +
  monitor("schematics/screen_network.png", 2.0, rows[3]) +
  # a protein
  protein(3.07, rows[1]) +
  # two drugs that act on a protein
  capsule(2.87, rows[2] + 0.22, angle = pi / 4, len = 0.12, r = 0.055, linewidth = 0.4) +
  capsule(3.27, rows[2] + 0.22, angle = -pi / 4, fill = "grey60", len = 0.12, r = 0.055, linewidth = 0.4) +
  fig_text(3.07, rows[2] + 0.22, "+", size = 10) +
  fig_arrow(2.91, rows[2] + 0.10, 3.0, rows[2] - 0.02, linewidth = 0.4, length = 0.04) +
  fig_arrow(3.23, rows[2] + 0.10, 3.14, rows[2] - 0.02, linewidth = 0.4, length = 0.04) +
  protein(3.07, rows[2] - 0.28, s = 0.16, linewidth = 0.3) +
  # a network of proteins, the middle row shows them as proteins
  annotate("segment", x = net$x[net_links$from], y = net$y[net_links$from],
           xend = net$x[net_links$to], yend = net$y[net_links$to], linewidth = 0.3) +
  geom_point(data = net[net$k != 2, ], aes(x, y), size = 1.3) +
  lapply(which(net$k == 2), function(i) protein(net$x[i], net$y[i], s = 0.09, linewidth = 0.2, halo = TRUE)) +
  coord_fixed(xlim = c(0.2, 3.45), ylim = c(0, 3.1), expand = FALSE) +
  theme_void()

# (b) the data and the model
# a cell with an irregular outline, a nucleus and a few organelles, s scales its size
cell <- function(x, y, phase = 0, s = 1, n = 120) {
  a <- seq(0, 2 * pi, length.out = n)
  radius <- 0.27 * s * (1 + 0.16 * cos(2 * a + phase) + 0.12 * cos(3 * a + 2 * phase) + 0.07 * cos(5 * a + 3 * phase))
  outline <- data.frame(x = x + radius * cos(a), y = y + radius * sin(a))
  nx <- x + 0.06 * s * cos(phase)
  ny <- y + 0.05 * s * sin(phase)
  nucleus <- data.frame(x = nx + 0.11 * s * cos(a), y = ny + 0.085 * s * sin(a + 0.4))
  dots <- data.frame(x = x + 0.17 * s * cos(phase + c(2.2, 3.1, 4.0, 4.9)),
                     y = y + 0.15 * s * sin(phase + c(2.2, 3.1, 4.0, 4.9)))
  list(geom_polygon(data = outline, aes(x, y), fill = "#e6eef8", colour = "grey20", linewidth = 0.5 * sqrt(s)),
       geom_polygon(data = nucleus, aes(x, y), fill = "#9fb3cf", colour = "grey20", linewidth = 0.35 * sqrt(s)),
       annotate("point", x = nx + 0.02 * s, y = ny + 0.01 * s, size = 0.7 * s, colour = "grey20"),
       geom_point(data = dots, aes(x, y), size = 0.35 * s, colour = "grey35"))
}
times <- data.frame(x = c(1.6, 2.1, 2.6), label = c("6 h", "24 h", "48 h"))
levels <- seq(2.0, 0.85, length.out = 6)
proteins <- expand.grid(k = 1:6, t = 1:3)
proteins$x <- times$x[proteins$t]
proteins$y <- levels[proteins$k]
# temporal dependencies between consecutive time points
links <- data.frame(t = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2),
                    from = c(1, 1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 4, 6, 6),
                    to = c(1, 3, 4, 2, 4, 6, 5, 2, 2, 5, 1, 4, 3, 6))
bars <- data.frame(x = 3.28 + 0.14 * 0:3, h = c(0.3, 0.48, 0.16, 0.38))

fig1_b <- ggplot() +
  cell(0.36, 1.83, phase = 0.3) + cell(0.66, 1.42, phase = 2.1) + cell(0.36, 1.0, phase = 4.0) +
  fig_text(0.5, 0.33, "Cell lines") +
  fig_arrow(1.08, 1.42, 1.4, 1.42) +
  # drug perturbations
  capsule(2.1, 2.82) +
  fig_text(2.38, 2.82, "Drug perturbations", hjust = 0) +
  fig_arrow(2.1, 2.6, 2.1, 2.27) +
  fig_text(2.19, 2.45, "alpha * ',' ~ beta", parse = TRUE, hjust = 0) +
  annotate("path", x = c(1.52, 1.52, 2.68, 2.68), y = c(2.15, 2.22, 2.22, 2.15), linewidth = 0.5) +
  # protein expression at the three time points
  annotate("segment", x = times$x[links$t], y = levels[links$from], xend = times$x[links$t + 1], yend = levels[links$to],
           colour = "grey55", linewidth = 0.3) +
  geom_point(data = proteins, aes(x, y, colour = factor(t)), size = 1.9) +
  scale_colour_manual(values = unname(fig_colors[c("blue", "orange", "green")]), guide = "none") +
  fig_text(times$x, 0.62, times$label, size = 8) +
  fig_text((times$x[-1] + times$x[-3]) / 2, 0.64, "gamma", parse = TRUE) +
  fig_text(2.1, 0.33, "Protein expression") +
  fig_arrow(2.84, 1.42, 3.2, 1.42) +
  # phenotype
  geom_rect(data = bars, aes(xmin = x, xmax = x + 0.1, ymin = 1.15, ymax = 1.15 + h), fill = "grey30") +
  annotate("segment", x = 3.24, xend = 3.84, y = 1.15, yend = 1.15, linewidth = 0.4) +
  fig_text(3.54, 0.33, "Phenotype") +
  coord_fixed(xlim = c(-0.35, 3.95), ylim = c(0, 3.1), expand = FALSE) +
  theme_void()

save_fig(list(fig1_a, fig1_b), 1, width = 7.5, height = 3.35, widths = c(3.25, 4.3), tags = c("(a)", "(b)"))

##### Figure 5: causal graph #####
nodes <- data.frame(name = c("D", "P6", "P24", "P48", "U"),
                    x = c(2, 0.8, 2, 3.2, 2), y = c(1, 0, 0, 0, -1),
                    label = c("italic(D)", "italic(P)^6", "italic(P)^24", "italic(P)^48", "italic(U)"),
                    hidden = c(FALSE, FALSE, FALSE, FALSE, TRUE))
edges <- data.frame(from = c("D", "D", "D", "P6", "P24", "U", "U", "U"),
                    to = c("P6", "P24", "P48", "P24", "P48", "P6", "P24", "P48"))
edges$hidden <- edges$from == "U"

# the arrows start and end at a distance r from the node centers
r <- 0.27
from <- nodes[match(edges$from, nodes$name), ]
to <- nodes[match(edges$to, nodes$name), ]
len <- sqrt((to$x - from$x)^2 + (to$y - from$y)^2)
edges$x <- from$x + r * (to$x - from$x) / len
edges$y <- from$y + r * (to$y - from$y) / len
edges$xend <- to$x - r * (to$x - from$x) / len
edges$yend <- to$y - r * (to$y - from$y) / len

u <- nodes[nodes$hidden, ]
angle <- seq(0, 2 * pi, length.out = 200)
circle <- data.frame(x = u$x + 0.22 * cos(angle), y = u$y + 0.22 * sin(angle))

# a dashed line would also dash the outline of its arrow head and could end in a cut off dash
# so each edge from U is drawn as a row of dashes of equal length, and the last dash carries a solid arrow head
head <- arrow(length = unit(0.07, "in"), type = "closed")
dashes <- do.call(rbind, lapply(which(edges$hidden), function(i) {
  e <- edges[i, ]
  len <- sqrt((e$xend - e$x)^2 + (e$yend - e$y)^2)
  dir <- c(e$xend - e$x, e$yend - e$y) / len
  # the dashes are 0.06 long with gaps of about 0.05, the last dash ends at the base of the arrow head
  start <- seq(0, len - 0.12, length.out = round((len - 0.01) / 0.11))
  data.frame(x = e$x + start * dir[1], y = e$y + start * dir[2],
             xend = e$x + (start + 0.06) * dir[1], yend = e$y + (start + 0.06) * dir[2], last = start == max(start))
}))
# the last dash is extended to the tip of the arrow
dashes$xend[dashes$last] <- edges$xend[edges$hidden]
dashes$yend[dashes$last] <- edges$yend[edges$hidden]

dag <- ggplot() +
  geom_segment(data = edges[!edges$hidden, ], aes(x, y, xend = xend, yend = yend), linewidth = 0.5, arrow = head) +
  geom_segment(data = dashes[!dashes$last, ], aes(x, y, xend = xend, yend = yend), colour = "grey50", linewidth = 0.5) +
  geom_segment(data = dashes[dashes$last, ], aes(x, y, xend = xend, yend = yend), colour = "grey50", linewidth = 0.5,
               arrow = head) +
  geom_path(data = circle, aes(x, y), colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_text(data = nodes, aes(x, y, label = label, colour = hidden), parse = TRUE,
            family = fig_family, size = 11 / .pt) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "grey50"), guide = "none") +
  coord_fixed(xlim = c(0.45, 3.55), ylim = c(-1.3, 1.25), expand = FALSE) +
  theme_void(base_size = 10, base_family = fig_family)

save_fig(dag, 5, width = 3.1, height = 2.55)


##### Striking images: two drawings without text in the style of Figure 1 #####
# all coordinates are in inches, both images are 7.5 in wide and 5 in high
time_colors <- unname(fig_colors[c("blue", "orange", "green")])

# a monitor with the screen (x0, x1, y0, y1), which shows two boxes in the style of the application
app_monitor <- function(x0, x1, y0, y1, bezel = 0.13) {
  cx <- (x0 + x1) / 2
  xs <- x0 + 0.66 * (x1 - x0)
  box <- function(a, b) list(
    geom_polygon(data = round_rect(a, b, y0 + 0.16, y1 - 0.3, 0.04), aes(x, y), fill = "white"),
    geom_polygon(data = round_rect(a, b, y1 - 0.5, y1 - 0.3, 0.04), aes(x, y), fill = "#3c8dbc"),
    annotate("rect", xmin = a, xmax = b, ymin = y1 - 0.5, ymax = y1 - 0.42, fill = "#3c8dbc"))
  # the lines of a list of results, the first three are highlighted
  rows <- data.frame(y = seq(y1 - 0.72, y0 + 0.38, length.out = 9),
                     len = c(0.8, 0.62, 0.9, 0.7, 0.55, 0.85, 0.6, 0.75, 0.5), top = 1:9 <= 3)
  list(
    # stand and frame
    annotate("rect", xmin = cx - 0.22, xmax = cx + 0.22, ymin = y0 - bezel - 0.5, ymax = y0, fill = "grey20"),
    geom_polygon(data = round_rect(cx - 0.85, cx + 0.85, y0 - bezel - 0.62, y0 - bezel - 0.46, 0.06), aes(x, y),
                 fill = "grey20"),
    geom_polygon(data = round_rect(x0 - bezel, x1 + bezel, y0 - bezel, y1 + bezel, 0.14), aes(x, y), fill = "grey20"),
    # page with a bar at the top
    annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = "#ecf0f5"),
    annotate("rect", xmin = x0, xmax = x1, ymin = y1 - 0.14, ymax = y1, fill = "#3c8dbc"),
    # a box for the network and a box with the list of results
    box(x0 + 0.16, xs - 0.08),
    box(xs + 0.08, x1 - 0.16),
    geom_rect(data = rows[rows$top, ], aes(xmin = xs + 0.14, xmax = x1 - 0.22, ymin = y - 0.07, ymax = y + 0.07),
              fill = "#cfe6f3"),
    geom_segment(data = rows, aes(x = xs + 0.2, xend = xs + 0.2 + len * (x1 - 0.5 - xs), y = y, yend = y),
                 colour = "grey60", linewidth = 0.9, lineend = "round")
  )
}

# the temporal network in the box (x0, x1, y0, y1), the proteins of the middle row are drawn as proteins
network <- function(x0, x1, y0, y1) {
  xs <- seq(x0, x1, length.out = 3)
  ys <- seq(y1, y0, length.out = 5)
  nodes <- expand.grid(k = 1:5, t = 1:3)
  nodes$x <- xs[nodes$t]
  nodes$y <- ys[nodes$k]
  links <- data.frame(t = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2),
                      from = c(1, 1, 2, 3, 4, 5, 5, 1, 2, 3, 3, 3, 4, 5),
                      to = c(1, 3, 3, 3, 2, 3, 5, 2, 1, 3, 1, 5, 4, 4))
  mid <- nodes[nodes$k == 3, ]
  list(annotate("segment", x = xs[links$t], y = ys[links$from], xend = xs[links$t + 1], yend = ys[links$to],
                colour = "grey55", linewidth = 0.55),
       geom_point(data = nodes[nodes$k != 3, ], aes(x, y, colour = factor(t)), size = 4),
       scale_colour_manual(values = time_colors, guide = "none"),
       lapply(1:3, function(i) {
         protein(mid$x[i], mid$y[i], s = 0.2, linewidth = 0.5, halo = TRUE, colour = time_colors[i])
       }))
}

# striking image 1: drugs, cells, and the application with the temporal network
# caption for the submission, one sentence per line
#   Drug-Prot quantifies how drugs and drug pairs change protein expression in breast cancer cell lines and how proteins depend on each other over time.
#   The image shows two drugs, a group of cells, and the application with the temporal protein network.
#   Blue, orange and green mark protein expression at 6, 24 and 48 hours after treatment.
striking1 <- ggplot() +
  app_monitor(3.2, 7.1, 1.45, 4.3) +
  network(3.75, 5.3, 1.95, 3.55) +
  capsule(0.95, 4.05, angle = pi / 4, len = 0.55, r = 0.24, linewidth = 1) +
  capsule(1.95, 4.05, angle = -pi / 4, fill = "grey60", len = 0.55, r = 0.24, linewidth = 1) +
  cell(0.95, 2.35, phase = 0.3, s = 2) + cell(1.8, 1.6, phase = 2.1, s = 2) + cell(0.9, 1.0, phase = 4.0, s = 2) +
  coord_fixed(xlim = c(0, 7.5), ylim = c(0, 5), expand = FALSE) +
  theme_void()
save_fig(striking1, "StrikingImage1", width = 7.5, height = 5)

# striking image 2: the model of Figure 1(b) without text
# caption for the submission, one sentence per line
#   Drug-Prot models the effect of drug perturbations on protein expression in cell lines at three time points and the dependence of each protein on the proteins at the preceding time point.
#   The image shows cell lines, a drug perturbation, the temporal protein network, and a phenotype such as the IC50.
#   Blue, orange and green mark protein expression at 6, 24 and 48 hours after treatment.
bars <- data.frame(x = 6.1 + 0.26 * 0:3, h = c(0.55, 0.9, 0.3, 0.7))
striking2 <- ggplot() +
  cell(0.85, 3.1, phase = 0.3, s = 1.9) + cell(1.45, 2.3, phase = 2.1, s = 1.9) + cell(0.85, 1.5, phase = 4.0, s = 1.9) +
  fig_arrow(2.15, 2.3, 2.7, 2.3, linewidth = 1, length = 0.13) +
  # drug perturbations
  capsule(4.0, 4.5, len = 0.42, r = 0.19, linewidth = 1) +
  fig_arrow(4.0, 4.12, 4.0, 3.72, linewidth = 1, length = 0.13) +
  annotate("path", x = c(2.9, 2.9, 5.1, 5.1), y = c(3.48, 3.6, 3.6, 3.48), linewidth = 0.9) +
  # protein expression at the three time points
  network(3.1, 4.9, 1.2, 3.25) +
  fig_arrow(5.35, 2.3, 5.9, 2.3, linewidth = 1, length = 0.13) +
  # phenotype
  geom_rect(data = bars, aes(xmin = x, xmax = x + 0.18, ymin = 1.85, ymax = 1.85 + h), fill = "grey30") +
  annotate("segment", x = 6.02, xend = 7.12, y = 1.85, yend = 1.85, linewidth = 0.7) +
  coord_fixed(xlim = c(0, 7.5), ylim = c(0, 5), expand = FALSE) +
  theme_void()
save_fig(striking2, "StrikingImage2", width = 7.5, height = 5)
