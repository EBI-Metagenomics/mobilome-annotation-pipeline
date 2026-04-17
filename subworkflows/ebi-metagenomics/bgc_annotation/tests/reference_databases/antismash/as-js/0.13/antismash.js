/*! antismash-js, version: 0.13.7 */
/*! For license information please see antismash.js.LICENSE.txt */
var viewer;
(() => {
    var t = {
            220: (t, e) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.clipboardCopyConstruct = e.copyToClipboard = void 0),
                    (e.copyToClipboard = function () {
                        const t = $(this).attr("data-seq") || "";
                        try {
                            navigator.clipboard.writeText(t);
                        } catch (t) {
                            return void alert(
                                "A browser permissions error occured while trying to copy the sequence to the clipboard",
                            );
                        }
                        alert("Sequence copied to clipboard");
                    }),
                    (e.clipboardCopyConstruct = function (t) {
                        return `<span class="clipboard-copy" data-seq="${t}">Copy to clipboard</span>`;
                    }));
            },
            225: (t, e) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.init = void 0));
                let n = null;
                function r(t) {
                    const e = $(`#${$(this).attr("id")}-tooltip`);
                    if ((n && (clearTimeout(e.data("timeout")), n.hide()), (n = e), "none" !== e.css("display")))
                        return void e.hide();
                    let r = $(this).offset();
                    void 0 === r && (r = { left: 0, top: 0 });
                    let i = setTimeout(() => e.slideUp("fast"), 5e3);
                    e.css("top", r.top + 10)
                        .css("left", r.left + 5)
                        .show()
                        .click(function () {
                            $(this).hide();
                        })
                        .data("timeout", i)
                        .mouseover(() => clearTimeout(e.data("timeout")))
                        .mouseout(() => {
                            (clearTimeout(e.data("timeout")),
                                (i = setTimeout(() => e.slideUp("fast"), 5e3)),
                                e.data("timeout", i));
                        });
                }
                e.init = function (t) {
                    $(`#${t} .clusterblast-orf`).each(function () {
                        const e = $(this);
                        (!(function (t, e) {
                            const n = t.attr("id"),
                                r = $("<div>");
                            (r.addClass("clusterblast-locustag"), r.attr("id", `${n}-label`));
                            const i = t.attr("locus_tag");
                            (i ? r.text(i) : r.text("unknown"),
                                $(`#${e}`).append(r),
                                t
                                    .mouseover((e) => {
                                        let i = t.offset();
                                        (void 0 === i && (i = { left: 0, top: 32 }),
                                            r.css("top", i.top - 32),
                                            r.css("left", i.left),
                                            $(`#${n}-label`).show());
                                    })
                                    .mouseout((t) => {
                                        $(`#${n}-label`).hide();
                                    }));
                        })(e, t),
                            (function (t, e) {
                                const n = t.attr("id"),
                                    i = $("<div>");
                                (i.addClass("clusterblast-tooltip"),
                                    i.attr("id", `${n}-tooltip`),
                                    i.html((t.attr("description") || "").replace("[br]", "<br>")),
                                    $(`#${e}`).append(i),
                                    t.click(r));
                            })(e, t));
                    });
                };
            },
            26: (t, e) => {
                "use strict";
                function n(t) {
                    (t.toggleClass("expanded"),
                        t.hasClass("expanded") ? t.next().css("display", "block") : t.next().removeAttr("style"));
                }
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.toggleCollapserHandler = e.toggleCollapser = void 0),
                    (e.toggleCollapser = n),
                    (e.toggleCollapserHandler = function (t) {
                        n($(this));
                    }));
            },
            322: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.setComparisonData = void 0));
                const r = n(86),
                    i = n(677),
                    a = n(3);
                let o = null,
                    s = null;
                function l(t, e, n, r, i) {
                    const a = n + r / 2,
                        o = n - r / 2,
                        s = n;
                    if (i) {
                        if (1 === t.strand) {
                            const n = Math.floor(e(t.start)),
                                i = Math.floor(Math.min(e(t.end) + r / 2, n));
                            return `${n},${a} ${i},${a} ${Math.floor(e(t.end))},${s} ${i},${o} ${n},${o}`;
                        }
                        if (-1 === t.strand) {
                            const n = Math.floor(e(t.start)),
                                i = Math.floor(e(t.end)),
                                l = Math.floor(Math.max(e(t.start) - r / 2, i));
                            return `${n},${s} ${l},${a} ${i},${a} ${i},${o} ${l},${o}`;
                        }
                        return `${t.start},${a} ${t.end},${a} ${t.end},${o} ${t.start},${o}`;
                    }
                    if (1 === t.strand) {
                        const n = Math.floor(e(t.start)),
                            i = Math.floor(Math.max(e(t.end) - r / 2, n));
                        return `${n},${a} ${i},${a} ${Math.floor(e(t.end))},${s} ${i},${o} ${n},${o}`;
                    }
                    if (-1 === t.strand) {
                        const n = Math.floor(e(t.start)),
                            i = Math.floor(e(t.end)),
                            l = Math.floor(Math.min(e(t.start) + r / 2, i));
                        return `${n},${s} ${l},${a} ${i},${a} ${i},${o} ${l},${o}`;
                    }
                    return `${t.start},${a} ${e(t.end)},${a} ${e(t.end)},${o} ${e(t.start)},${o}`;
                }
                function c(t, e) {
                    const n = window.location.hash.substring(1),
                        c = document.getElementById(`comparison-${e}-${n}-selector`);
                    if (!c || !c.value) return;
                    const u = $(`#${t} * .heat-row-${e}`);
                    (u.off("click").click(function () {
                        ($(`.heat-row-${e}`).css("border", ""), $(this).css("border", "1px solid black"));
                        const t = $(this),
                            u = `comparison-${e}-${n}-svg`;
                        $(`#${u}`).remove();
                        const h = (0, i.select)(`#comparison-${e}-${n}`)
                                .append("svg")
                                .attr("width", "100%")
                                .attr("id", u),
                            d = o[e][c.value],
                            g = t.attr("data-accession");
                        g &&
                            ((function (t, e, n, i) {
                                if (null === i || null === s) return;
                                const o = 54.5,
                                    c = 109.5,
                                    u = s.end - s.start,
                                    h = s.idx,
                                    d = 800,
                                    g = (0, r.scaleLinear)().domain([s.start, s.end]).range([0, d]);
                                let f = (0, r.scaleLinear)().domain([i.start, i.end]).range([0, d]);
                                (i.reverse && (f = (0, r.scaleLinear)().domain([i.end, i.start]).range([0, d])),
                                    t.attr("viewport", `0 164 0 ${u}`).attr("height", 164).attr("width", d),
                                    t
                                        .append("text")
                                        .attr("x", 0)
                                        .attr("y", 20)
                                        .attr("text-anchor", "start")
                                        .text("Query"),
                                    t
                                        .append("text")
                                        .attr("x", 0)
                                        .attr("y", 155)
                                        .attr("text-anchor", "start")
                                        .text(`Reference: ${e.indexOf("-mibig-") > -1 ? n.split(":")[0] : n}`),
                                    t
                                        .append("line")
                                        .attr("x1", 0)
                                        .attr("y1", o)
                                        .attr("x2", u)
                                        .attr("y2", o)
                                        .attr("class", "centerline"),
                                    t
                                        .append("line")
                                        .attr("x1", 0)
                                        .attr("y1", c)
                                        .attr("x2", u)
                                        .attr("y2", c)
                                        .attr("class", "centerline"),
                                    t
                                        .selectAll("line.link")
                                        .data(i.links)
                                        .enter()
                                        .append("line")
                                        .attr("x1", (t) => g(t.query_loc))
                                        .attr("y1", 59.5)
                                        .attr("x2", (t) => f(t.subject_loc))
                                        .attr("y2", 104.5)
                                        .attr("class", "link")
                                        .attr("class", "centerline"));
                                const p = t
                                    .selectAll("g.cc-svg-orf-group")
                                    .data(s.orfs)
                                    .enter()
                                    .append("g")
                                    .attr("class", "cc-svg-orf-group");
                                (p
                                    .append("polygon")
                                    .attr("points", (t) => l(t, g, o, 10))
                                    .attr("class", "cc-svg-orf-bg")
                                    .attr("id", (t) => `${e}-${(0, a.tag_to_id)(t.locus_tag)}-bg`)
                                    .style("fill", "white"),
                                    p
                                        .append("polygon")
                                        .attr("points", (t) => l(t, g, o, 10))
                                        .attr("class", (t) => `cc-svg-orf svgene-type-${t.type}`)
                                        .attr("id", (t) => `${e}-${(0, a.tag_to_id)(t.locus_tag)}`)
                                        .attr("opacity", "1")
                                        .attr("data-locus", (t) => t.locus_tag),
                                    t
                                        .selectAll("text.svgene-locustag")
                                        .data(s.orfs)
                                        .enter()
                                        .append("text")
                                        .attr("x", (t) => (g(t.start) < 400 ? g(t.start) : g(t.end)))
                                        .attr("text-anchor", (t) => (g(t.start) < 400 ? "start" : "end"))
                                        .attr("y", 36.5)
                                        .attr("class", "cc-svg-locustag")
                                        .attr("id", (t) => `${e}-${(0, a.tag_to_id)(t.locus_tag)}-tag`)
                                        .text((t) => t.locus_tag));
                                const m = t
                                    .selectAll("g.cc-svg-reforf-group")
                                    .data(i.genes)
                                    .enter()
                                    .append("g")
                                    .attr("class", "cc-svg-reforf-group");
                                (m
                                    .append("polygon")
                                    .attr("points", (t) => l(t, f, c, 10, i.reverse))
                                    .attr("class", "cc-svg-reforf-bg")
                                    .attr("id", (t) => `${e}-ref${(0, a.tag_to_id)(t.locus_tag)}-bg`)
                                    .style("fill", "white"),
                                    m
                                        .append("polygon")
                                        .attr("points", (t) => l(t, f, c, 10, i.reverse))
                                        .attr("class", (t) => `cc-svg-reforf svgene-type-${t.function}`)
                                        .attr("id", (t) => `${e}-ref${(0, a.tag_to_id)(t.locus_tag)}`)
                                        .attr("data-locus", (t) => t.locus_tag)
                                        .style("opacity", (t) => (t.linked[`${h}`] ? "1" : "0.5")),
                                    t
                                        .selectAll("text.svgene-locustag")
                                        .data(i.genes)
                                        .enter()
                                        .append("text")
                                        .attr("x", (t) => (f(t.start) < 400 ? f(t.start) : f(t.end)))
                                        .attr("text-anchor", (t) => (f(t.start) < 400 ? "start" : "end"))
                                        .attr("y", 127.5)
                                        .attr("class", "cc-svg-locustag")
                                        .attr("id", (t) => `${e}-ref${(0, a.tag_to_id)(t.locus_tag)}-tag`)
                                        .text((t) => t.locus_tag));
                            })(h, u, g, d.reference_clusters[g]),
                            $(".cc-svg-orf")
                                .mouseover(function (t) {
                                    const e = $(this).attr("id");
                                    void 0 !== e && $(`#${e}-tag`).show();
                                })
                                .mouseout(function (t) {
                                    const e = $(this).attr("id");
                                    void 0 !== e && $(`#${e}-tag`).hide();
                                })
                                .click(function (t) {
                                    const e = $(this).attr("data-locus");
                                    e && (0, a.selectOrfsByLoci)([e], t.ctrlKey || t.metaKey);
                                }),
                            $(".cc-svg-reforf")
                                .mouseover(function (t) {
                                    const e = $(this).attr("id");
                                    void 0 !== e && $(`#${e}-tag`).show();
                                })
                                .mouseout(function (t) {
                                    const e = $(this).attr("id");
                                    void 0 !== e && $(`#${e}-tag`).hide();
                                }));
                    }),
                        u.first().click());
                }
                e.setComparisonData = function (t, e, n) {
                    ((o = e), (s = n));
                    for (const n in e) e.hasOwnProperty(n) && c(t, n);
                };
            },
            128: (t, e) => {
                "use strict";
                function n(t, e) {
                    !(function (t, e) {
                        ($(`#${e} * .body-details-section`).hide(),
                            $(`#${e} * .body-details-header-active`).toggleClass("body-details-header-active"),
                            t.addClass("body-details-header-active"),
                            $(`#${e} * .body-details-section.${t.attr("data-name")}`).show());
                    })($(this), t);
                    const n = $(`#${t} * .sidepanel-details-header.${$(this).attr("data-name")}`);
                    n.length > 0 && i(n, t);
                }
                function r(t, e) {
                    (i($(this), t), $(`#${t} * .body-details-header.${$(this).attr("data-name")}`).click());
                }
                function i(t, e) {
                    ($(`#${e} * .sidepanel-details-section`).hide(),
                        $(`#${e} * .sidepanel-details-header-active`).toggleClass("sidepanel-details-header-active"),
                        t.addClass("sidepanel-details-header-active"),
                        $(`#${e} * .sidepanel-details-section.${t.attr("data-name")}`).show());
                }
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.setupDetails = void 0),
                    (e.setupDetails = function (t) {
                        ($(".sidepanel-details-header").off("click"), $(".sidepanel-details-section").hide());
                        for (const e of t)
                            $(`#${e} * .sidepanel-details-header`)
                                .click($.proxy(r, null, e))
                                .first()
                                .trigger("click");
                        ($(".body-details-header").off("click"), $(".body-details-section").hide());
                        for (const e of t)
                            $(`#${e} * .body-details-header`)
                                .click($.proxy(n, null, e))
                                .first()
                                .trigger("click");
                    }));
            },
            501: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.actualDrawDomainBubbleData = e.drawDomainBubbleData = void 0));
                const r = n(677),
                    i = n(3),
                    a = 15,
                    o = -a,
                    s = m(0, o, 2 * a) - a,
                    l = 6,
                    c = a / 3,
                    u = 2 * a,
                    h = 16,
                    d = 10,
                    g = 20,
                    f = 20;
                let p = null;
                function m(t, e, n) {
                    return 0 === t
                        ? Math.round(Math.sqrt(n * n - e * e))
                        : 0 === e
                          ? Math.round(Math.sqrt(n * n - t * t))
                          : Math.round(Math.sqrt(t * t + e * e));
                }
                function v(t, e) {
                    const n = e;
                    let r = n.y - 1.5 * l,
                        i = n.y + 1.5 * l,
                        a = n.x - l,
                        o = a + 2 * l;
                    const s = l - 1;
                    return (
                        t
                            .append("polygon")
                            .attr("class", "bubble-domain-terminal-docking")
                            .attr(
                                "points",
                                `${a},${r} ${a + s},${n.y} ${a},${i} ${a + s},${i} ${o},${n.y} ${a + s},${r} ${a},${r}`,
                            ),
                        "end" === e.terminalDocking ? ((a = n.x - l), (o = n.x)) : ((a = n.x), (o = n.x + l)),
                        (r = n.y - l / 3),
                        (i = n.y + l / 3),
                        t
                            .append("polygon")
                            .attr("class", "bubble-domain-terminal-docking")
                            .attr("points", `${a},${r} ${o},${r} ${o},${i} ${a},${i} ${a},${r}`),
                        t
                    );
                }
                function y(t, e, n) {
                    if (!e.iterative)
                        return (
                            t
                                .append("line")
                                .attr("x1", e.start)
                                .attr("x2", e.end)
                                .attr("y1", n)
                                .attr("y2", n)
                                .attr("class", "bubble-module-line" + (e.nonElongating ? "-non-elongating" : "")),
                            t
                        );
                    const r = n + 2.5 * a,
                        i = n + (r - n) / 2,
                        o = (e.end, r - i + 2),
                        s = e.start + 1.5 * a,
                        l = e.end - 1.5 * a,
                        c = e.start + 0.75 * (e.end - e.start);
                    return (
                        t
                            .append("path")
                            .attr("class", "bubble-module-line")
                            .attr(
                                "d",
                                `M ${
                                    e.start
                                },${i}\n                  A ${o} ${o} 0 0 1 ${s},${n}\n                  L ${l},${n}\n                  A ${o} ${o} 0 0 1 ${
                                    e.end
                                },${i}\n                  A ${o} ${o} 0 0 1 ${l},${r}\n                  L ${c},${r}\n                  L ${
                                    c + 10
                                },${r - 10}\n                  M ${
                                    c - 10
                                },${r}\n                  L ${s},${r}\n                  A ${o} ${o} 0 0 1 ${
                                    e.start
                                },${i}`,
                            )
                            .attr("fill", "none"),
                        t
                    );
                }
                function b(t, e) {
                    if (!document.getElementById(`${t}-domain-bubble-svg-container`)) return;
                    const n = `${t}-domain-bubble-svg`;
                    if (($(`#${n}`).remove(), !e)) return;
                    const h = [],
                        p = e.modules,
                        b = (function (t) {
                            const e = [];
                            for (const n of t) for (const t of n.domains) e.push({ dom: t, mod: t.modifier });
                            e.push({ dom: null, mod: !1 });
                            let n = 0,
                                r = 0,
                                i = 0;
                            for (const t of e) {
                                if ((t.dom && (t.dom.level = 0), t.mod)) {
                                    ((r += 1), (n += 1));
                                    continue;
                                }
                                if (0 === r) {
                                    n += 1;
                                    continue;
                                }
                                const a = Math.floor(r / 2);
                                a + (r % 2) > i && (i = a + (r % 2));
                                for (let t = 0; t < a; t++)
                                    ((e[n - r + t].dom.level = t + 1), (e[n - 1 - t].dom.level = t + 1));
                                (r % 2 == 1 && (e[n - r + a].dom.level = a + 1), (r = 0), (n += 1));
                            }
                            return i + 1;
                        })(p),
                        A = a * (b + 1) + 2 * g,
                        _ = (0, r.select)(`#${t}-domain-bubble-svg-container`)
                            .append("svg")
                            .attr("width", "100%")
                            .attr("height", A + 4.5 * a)
                            .attr("id", n)
                            .attr("class", "domain-bubble-container");
                    let C = a,
                        S = A + o * p[0].domains[0].level,
                        M = C,
                        T = 0,
                        k = 0,
                        R = null,
                        N = "",
                        L = 0;
                    for (const t of p) {
                        (T > 0 && (C += c), (t.number = ++T), t.complete && (t.completeNumber = ++k));
                        let e = 0;
                        for (const n of t.domains)
                            (R &&
                                (R.cds !== n.cds
                                    ? (0 === e && (C -= c),
                                      h.push({ end: C + R.radius, name: N, scale: 1, start: M }),
                                      (C += u),
                                      (M = C),
                                      (L = 0))
                                    : (L =
                                          R && R.special && n.modifier
                                              ? a - m(0, o, l + a)
                                              : R.level !== n.level
                                                ? s
                                                : R.radius),
                                R.cds === n.cds && (n.terminalDocking || R.terminalDocking) && (L -= c)),
                                (C += L),
                                n.special ? (n.radius = l) : (n.radius = a),
                                (C += n.radius),
                                (n.x = C),
                                R && R.level !== n.level && (R.level > n.level ? (S -= o) : (S += o)),
                                n.terminalDocking && (S = A),
                                (n.y = S),
                                (R = n),
                                (N = n.cds),
                                (e += 1));
                        const n = t.domains[0],
                            r = t.domains[t.domains.length - 1];
                        ((t.start = n.x - n.radius), (t.end = r.x + r.radius), (L = 0));
                    }
                    (R && (C += R.radius), h.push({ end: C, name: N, scale: 1, start: M }));
                    const P = (0, r.select)(".bubble-tooltip");
                    function D(t) {
                        const e = (0, r.select)(this);
                        if (t.terminalDocking) return void v(e, t);
                        (e
                            .append("circle")
                            .attr("class", (t) => `${t.css}`)
                            .attr("cx", (t) => t.x)
                            .attr("cy", (t) => t.y)
                            .attr("r", (t) => (t.special ? l : a)),
                            e
                                .append("text")
                                .attr("x", (t) => t.x)
                                .attr("y", (t) => t.y + a / 3)
                                .attr("text-anchor", "middle")
                                .text((t) => (t.special ? "" : t.name)));
                        const n = Math.sqrt((a * a) / 2),
                            i = e
                                .filter((t) => t.inactive)
                                .append("g")
                                .attr("class", "bubble-domain-inactive");
                        (i
                            .append("line")
                            .attr("x1", (t) => t.x - n)
                            .attr("x2", (t) => t.x + n)
                            .attr("y1", (t) => t.y - n)
                            .attr("y2", (t) => t.y + n),
                            i
                                .append("line")
                                .attr("x1", (t) => t.x - n)
                                .attr("x2", (t) => t.x + n)
                                .attr("y1", (t) => t.y + n)
                                .attr("y2", (t) => t.y - n));
                    }
                    T = 0;
                    for (const t of p) {
                        const e = _.append("g")
                            .attr("class", `module-group-${t.number}`)
                            .selectAll("g.domain-group")
                            .data(t.domains)
                            .enter()
                            .append("g")
                            .attr("class", "g.domain-group");
                        (e.each(D),
                            t.complete || e.style("opacity", "50%"),
                            e
                                .on("click", (t) =>
                                    P.text(t.description)
                                        .style("display", "block")
                                        .style("top", `${r.event.pageY + 20}px`)
                                        .style("left", `${r.event.pageX + 20}px`),
                                )
                                .on("mouseout", () => P.style("display", "none")));
                    }
                    S = A;
                    const I = _.append("g")
                        .attr("class", "module-labels")
                        .selectAll("g.bubble-module-label")
                        .data(p.filter((t) => t.complete))
                        .enter()
                        .append("g")
                        .attr("class", "bubble-module-label");
                    (I.each(function (t) {
                        return y((0, r.select)(this), t, S + 1.5 * a);
                    })
                        .append("text")
                        .attr("x", (t) => (t.start + t.end) / 2)
                        .attr("y", S + 2.5 * a)
                        .attr("class", "bubble-module-number")
                        .attr("text-anchor", "middle")
                        .text((t) => `M ${t.completeNumber}`),
                        I.append("text")
                            .attr("x", (t) => (t.start + t.end) / 2)
                            .attr("y", S + 3.5 * a)
                            .attr("class", "bubble-module-monomer")
                            .attr("text-anchor", "middle")
                            .text((t) => t.polymer),
                        (S -= a * (b + 1)));
                    const E = _.append("g")
                        .attr("class", "bubble-genes")
                        .selectAll("line.bubble-gene-line")
                        .data(h)
                        .enter()
                        .append("g");
                    (E.attr("data-locus", (t) => t.name).attr("class", "bubble-gene"),
                        E.append("polygon")
                            .attr("class", "bubble-gene-arrow")
                            .attr("points", (t) => w(t.start, t.end, S)),
                        E.append("text")
                            .text((t) => t.name)
                            .attr("x", (t) => (t.start + t.end) / 2)
                            .attr("y", S)
                            .attr("class", "serif bubble-gene-label")
                            .style("font-size", "1px")
                            .each(x)
                            .style("font-size", (t) => (t.scale >= d ? `${t.scale}px` : "10px"))
                            .style("fill", (t) => (t.scale >= d ? "white" : "black"))
                            .attr("text-anchor", (t, e) => {
                                if (t.scale < d) {
                                    if (0 === e) return "start";
                                    if (e === h.length - 1) return "end";
                                }
                                return "middle";
                            })
                            .attr("x", (t, e) => {
                                const n = Math.floor((t.start + t.end) / 2);
                                return t.scale < d ? (0 === e ? t.start : e === h.length - 1 ? t.end : n) : n - f / 5;
                            })
                            .attr("y", (t, e) =>
                                t.scale < d ? S - 1.5 * g + ((e + 1) % 2) * 15 : Math.floor(S + t.scale / 3),
                            ),
                        _.attr("width", C + 2 * a),
                        $(".bubble-gene")
                            .off("click")
                            .click(function (t) {
                                $(`#${(0, i.locusToFullId)($(this).attr("data-locus") || "none")}-svgeneorf`).trigger(
                                    t,
                                );
                            }));
                }
                function x(t) {
                    if (!this || !this.parentNode) return void (t.scale = h);
                    const e = this.getBBox(),
                        n = this.parentNode.getBBox();
                    return (n.height, (t.scale = Math.min(h, Math.round((n.width - f) / e.width))), t.scale);
                }
                function w(t, e, n) {
                    const r = n - g / 2,
                        i = r + g;
                    let a = f;
                    0.2 * (e - t) < a && (a = Math.floor(0.2 * (e - t)));
                    const o = Math.floor((i - r) / 5),
                        s = e - a;
                    return `${t},${r}\n            ${s},${r}\n            ${s},${
                        r - o
                    }\n            ${e},${n}\n            ${s},${
                        i + o
                    }\n            ${s},${i}\n            ${t},${i}\n            ${t},${r}`;
                }
                ((e.drawDomainBubbleData = function (t, e) {
                    ((p = e),
                        $(`#${t}-domain-bubble-select`)
                            .off("change")
                            .change(function () {
                                const e = "" + $(this).val();
                                if (p && e && p[e]) {
                                    const n = {};
                                    (b(t, p[e]),
                                        (n.showNonElongating = p[e].modules.some((t) => t.nonElongating)),
                                        (function (t, e) {
                                            const n = "bubble-legend";
                                            if ($(`#${n}`).length)
                                                return void $(`#${n}`).detach().appendTo(`#${t}-bubble-legend`);
                                            const i = 2 * a + 4,
                                                o = 2 * a + 4,
                                                s = a + 2,
                                                c = (0, r.select)(`#${t}-bubble-legend`).append("div").attr("id", n);
                                            (c
                                                .append("div")
                                                .attr("class", "bubble-legend-icon")
                                                .append("svg")
                                                .attr("width", o)
                                                .attr("height", i)
                                                .append("circle")
                                                .attr("class", "jsdomain-other module-bubble")
                                                .attr("cx", s)
                                                .attr("cy", s)
                                                .attr("r", a),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("Domain in a complete module"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon")
                                                    .append("svg")
                                                    .attr("width", o)
                                                    .attr("height", i)
                                                    .append("circle")
                                                    .attr("class", "jsdomain-other module-bubble")
                                                    .attr("cx", s)
                                                    .attr("cy", s)
                                                    .attr("r", l),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("Special domain (e.g. trans-AT docking domains)"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon")
                                                    .append("svg")
                                                    .attr("width", o)
                                                    .attr("height", i)
                                                    .append("circle")
                                                    .attr("class", "jsdomain-other module-bubble")
                                                    .attr("cx", s)
                                                    .attr("cy", s)
                                                    .attr("r", a)
                                                    .style("opacity", "50%"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("Domain in an incomplete module or outside modules"),
                                                v(
                                                    c
                                                        .append("div")
                                                        .attr("class", "bubble-legend-icon")
                                                        .append("svg")
                                                        .attr("width", o)
                                                        .attr("height", i),
                                                    { x: o / 2, y: i / 2, terminalDocking: "start" },
                                                ).style("opacity", "50%"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("N-terminal docking domains"),
                                                v(
                                                    c
                                                        .append("div")
                                                        .attr("class", "bubble-legend-icon")
                                                        .append("svg")
                                                        .attr("width", o)
                                                        .attr("height", i),
                                                    { x: o / 2, y: i / 2, terminalDocking: "end" },
                                                ).style("opacity", "50%"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("C-terminal docking domains"));
                                            const u = Math.floor(Math.sqrt((a * a) / 2)),
                                                h = s + u,
                                                d = s - u,
                                                g = c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon")
                                                    .append("svg")
                                                    .attr("width", o)
                                                    .attr("height", i)
                                                    .append("g")
                                                    .attr("class", "bubble-domain-inactive");
                                            (g.append("line").attr("x1", d).attr("y1", d).attr("x2", h).attr("y2", h),
                                                g
                                                    .append("line")
                                                    .attr("x1", d)
                                                    .attr("y1", h)
                                                    .attr("x2", h)
                                                    .attr("y2", d),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("Marks domains that are predicted to be inactive"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon")
                                                    .append("svg")
                                                    .attr("width", 4 * a)
                                                    .attr("height", i)
                                                    .append("polygon")
                                                    .attr("class", "bubble-gene-arrow")
                                                    .attr("points", w(a / 3, 4 * a - a / 3, s)),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("The gene/CDS feature containing the domains"));
                                            const f = c
                                                .append("div")
                                                .attr("class", "bubble-legend-icon")
                                                .append("svg")
                                                .attr("width", 3 * a)
                                                .attr("height", i)
                                                .append("g")
                                                .attr("class", "bubble-module-label");
                                            if (
                                                (f
                                                    .append("line")
                                                    .attr("x1", 2)
                                                    .attr("x2", 3 * a - 2)
                                                    .attr("y1", 4)
                                                    .attr("y2", 4)
                                                    .attr("class", "bubble-module-line"),
                                                f
                                                    .append("text")
                                                    .attr("x", 1.75 * a)
                                                    .attr("y", 1.75 * a)
                                                    .attr("class", "bubble-module-text")
                                                    .attr("text-anchor", "middle")
                                                    .text("M #"),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text(
                                                        "The extent of a module, with the module number and the predicted monomer for that module",
                                                    ),
                                                e.showNonElongating)
                                            ) {
                                                const t = c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon bubble-legend-non-elongating")
                                                    .append("svg")
                                                    .attr("width", 3 * a)
                                                    .attr("height", i)
                                                    .append("g")
                                                    .attr("class", "bubble-module-label");
                                                (t
                                                    .append("line")
                                                    .attr("x1", 2)
                                                    .attr("x2", 3 * a - 2)
                                                    .attr("y1", 4)
                                                    .attr("y2", 4)
                                                    .attr("class", "bubble-module-line-non-elongating"),
                                                    t
                                                        .append("text")
                                                        .attr("x", 1.75 * a)
                                                        .attr("y", 1.75 * a)
                                                        .attr("class", "bubble-module-text")
                                                        .attr("text-anchor", "middle")
                                                        .text("M #"),
                                                    c
                                                        .append("div")
                                                        .attr(
                                                            "class",
                                                            "bubble-legend-text bubble-legend-non-elongating",
                                                        )
                                                        .text(
                                                            "The extent of a module predicted to be non-elongating, though modification domains may still apply",
                                                        ));
                                            }
                                            (y(
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-icon")
                                                    .append("svg")
                                                    .attr("width", 6 * a)
                                                    .attr("height", 2.75 * a)
                                                    .attr("viewbox", `0 0 ${6 * a} ${4 * a}`)
                                                    .attr("transform", `scale(${i / (3.5 * a)} ${i / (3.5 * a)})`)
                                                    .append("g")
                                                    .attr("class", "bubble-module-label"),
                                                {
                                                    complete: !0,
                                                    completeNumber: 1,
                                                    domains: [],
                                                    end: 6 * a - 2,
                                                    iterative: !0,
                                                    number: 1,
                                                    polymer: "",
                                                    start: 2,
                                                },
                                                1,
                                            ),
                                                c
                                                    .append("div")
                                                    .attr("class", "bubble-legend-text")
                                                    .text("Module extents when predicted to be iterative"),
                                                $(".bubble-legend-icon").css("height", `${i}px`));
                                        })(t, n));
                                }
                            }),
                        $(`.body-details-header.${t}-nrps-pks-bubbles`)
                            .off(".firstClick")
                            .one("click.firstClick", () => {
                                $(`#${t}-domain-bubble-select`).change();
                            }),
                        0 === $(".bubble-tooltip").length &&
                            (0, r.select)("body").append("div").attr("class", "tooltip bubble-tooltip").text(""));
                }),
                    (e.actualDrawDomainBubbleData = b));
            },
            415: (t, e) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.downloadSvg = e.initDownloadButtons = void 0));
                const n = "http://www.w3.org/2000/svg",
                    r = "http://www.w3.org/1999/xlink",
                    i = "http://www.w3.org/2000/xmlns/";
                function a(t) {
                    if (!t.target) return !0;
                    const n = t.target,
                        r = n.getAttribute("data-tag"),
                        i = n.getAttribute("data-filename");
                    return !(!r || !i || ((0, e.downloadSvg)(r, i), 1));
                }
                function o(t, e) {
                    const n = getComputedStyle(t);
                    let r = "";
                    for (let t = 0, i = n.length; t < i; t++) {
                        const i = n[t],
                            a = n.getPropertyValue(i);
                        "cursor" !== i && a !== e.get(i) && (r += `${i}:${a};`);
                    }
                    t.setAttribute("style", r);
                }
                function s(t, e) {
                    if (!(null == t ? void 0 : t.hasChildNodes())) return;
                    let n = t.firstChild;
                    for (; n; )
                        (n.nodeType === Node.ELEMENT_NODE && "SCRIPT" !== n.nodeName && (e.push(n), s(n, e)),
                            (n = n.nextSibling));
                }
                ((e.initDownloadButtons = () => {
                    document.querySelectorAll(".download-svg").forEach((t) => {
                        (t.removeEventListener("click", a, !1), t.addEventListener("click", a, !1));
                    });
                }),
                    (e.downloadSvg = (t, e) => {
                        const a = document.querySelector(`#${t}`);
                        if (!a) return;
                        let l;
                        if (((l = "svg" === a.nodeName ? a : a.children.item(0)), !l || "svg" !== l.nodeName)) return;
                        const c = l.cloneNode(!0);
                        !(function (t) {
                            (t.hasAttribute("version") || t.setAttribute("version", "1.1"),
                                t.hasAttributeNS(i, "xmlns") || t.setAttributeNS(i, "xmlns", n),
                                t.hasAttributeNS(i, "xmlns:xlink") || t.setAttributeNS(i, "xmlns:xlink", r),
                                document.body.appendChild(t),
                                (function (t) {
                                    const e = document.createElementNS(n, "svg");
                                    document.body.appendChild(e);
                                    const r = new Map(),
                                        i = getComputedStyle(e);
                                    (Object.values(i).forEach((t) => {
                                        r.set(t, i.getPropertyValue(t));
                                    }),
                                        document.body.removeChild(e));
                                    const a = (function (t) {
                                        const e = [];
                                        return (e.push(t), s(t, e), e);
                                    })(t);
                                    let l = a.length;
                                    for (; l--; ) o(a[l], r);
                                    (t.removeAttribute("width"),
                                        t.removeAttribute("height"),
                                        t.style.removeProperty("width"),
                                        t.style.removeProperty("height"),
                                        t.style.removeProperty("inline-size"));
                                    const c = t.getAttribute("viewbox");
                                    c && c.split(" ").length > 3 && t.setAttribute("width", c.split(" ")[2]);
                                })(t),
                                document.body.removeChild(t));
                        })(c);
                        const u = new XMLSerializer().serializeToString(c),
                            h = `data:image/svg+xml;base64,${btoa(u)}`,
                            d = document.createElement("a");
                        ((d.href = h), (d.download = e));
                        const g = (t) => {
                            setTimeout(() => {
                                var e;
                                null === (e = null == t ? void 0 : t.target) ||
                                    void 0 === e ||
                                    e.removeEventListener("click", g);
                            }, 200);
                        };
                        (d.addEventListener("click", g, !1), d.click());
                    }));
            },
            457: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.initGeneTableHandler = void 0));
                const r = n(220),
                    i = n(3),
                    a = n(130),
                    o = "gt";
                function s(t, e, n) {
                    if (!n)
                        for (const e of t.querySelectorAll(`.${o}-matching-cell`))
                            e.classList.remove(`.${o}-matching-cell`);
                    (e % 2 == 1
                        ? (t.classList.remove("row-even"), t.classList.add("row-odd"))
                        : (t.classList.remove("row-odd"), t.classList.add("row-even")),
                        (t.style.display = ""));
                    const r = t.querySelectorAll(`td.${o}-cell-match`)[0];
                    r && n.length > 0 && ((r.style.display = n.length > 0 ? "" : "none"), (r.innerHTML = n));
                }
                e.initGeneTableHandler = function (t, e) {
                    const n = `#${o}-container-${t.anchor}`,
                        l = `#${o}-${t.anchor}`,
                        c = $(n),
                        u = c.find(`.${o}-scroll-container`),
                        h = $(l),
                        d = h.find("thead"),
                        g = c.find(`.${o}-no-hits`),
                        f = $(`#${o}-${t.anchor}-toggle`);
                    var p, m;
                    if (
                        ((p = f),
                        (m = $(`.${o}-zoom-toggle`)).off("change"),
                        p.on("change", () => {
                            const t = p.prop("checked");
                            (m.not(p).prop("checked", t), t ? (0, i.zoom_to_selection)() : (0, i.resetZoom)());
                        }),
                        0 === h.length)
                    )
                        return;
                    !(function (t, e, n) {
                        const r = $(`${t} tbody`);
                        if ($(r).children().length > 0) return;
                        let i = 1;
                        for (const t of e.orfs.sort((t, e) => t.start - e.start)) {
                            const e = (0, a.replaceWildcardsInText)(n.blast_template, t),
                                o = [
                                    `<td class="serif search-cell" data-type="name">${(0, a.replaceWildcardsInText)(
                                        n.selected_template,
                                        t,
                                    )}${t.locus_tag}</td>`,
                                    `<td class="search-cell" data-type="product">${t.product}</td>`,
                                    `<td class="gt-cell-numeric">${t.dna.length}</td>`,
                                    `<td class="gt-cell-numeric">${t.translation.length}</td>`,
                                    `<td class="search-cell" data-type="function"><span class="colour-square legend-type-${t.type}"></span>${t.type}</td>`,
                                    `<td class="seq-copy seq-copy-nt gt-cell-interactable" data-locus="${t.locus_tag}"><span class="button-like">Copy</span></td>`,
                                    `<td class="seq-copy seq-copy-aa gt-cell-interactable" data-locus="${t.locus_tag}"><span class="button-like">Copy</span></td>`,
                                    `<td class="gt-cell-link gt-cell-interactable">${e}</td>`,
                                    '<td class="gt-column-match gt-cell-match" style="display:none"></td>',
                                ];
                            (r.append(
                                `<tr class="gt-row row-${i % 2 ? "odd" : "even"}" data-locus="${t.locus_tag}">${o.join(
                                    "",
                                )}</tr>`,
                            ),
                                i++);
                        }
                    })(l, t, e);
                    let v = !0;
                    const y = new MutationObserver((t, e) => {
                        if (!v || 0 === u.length) return;
                        const n = $(h.find(".cds-selected-marker.active").first().closest("tr"));
                        if (!n.length) return;
                        u.scrollTop(0);
                        const r = n.position(),
                            i = d.position(),
                            a = d.height();
                        r && i && a && u.animate({ scrollTop: r.top - i.top - a }, "250");
                    });
                    h.find(".cds-selected-marker").each(function () {
                        y.observe(this, { attributeFilter: ["class"], attributes: !0 });
                    });
                    const b = $(`${l}-search`);
                    (b
                        .off(`input.${o}-search`)
                        .on("input", function (t) {
                            if ((t.stopPropagation(), !(this instanceof HTMLInputElement))) return;
                            const n = this.value;
                            $(`.${o}-matching-cell`).removeClass(`${o}-matching-cell`);
                            const r = h.find(`.${o}-column-match`);
                            if (0 === n.length) {
                                let t = 1;
                                (h.show(),
                                    h.find(`.${o}-row`).each(function () {
                                        s(this, t++, "");
                                    }),
                                    h.find(`.${o}-column-match`).hide(),
                                    g.hide(),
                                    f.prop("checked") && (0, i.resetZoom)());
                            } else
                                (!(function (t, e, n) {
                                    let r = 0;
                                    const a = [];
                                    let l = !1;
                                    for (const i of e.getElementsByTagName("tbody")[0].querySelectorAll("tr")) {
                                        const e = i.getAttribute("data-locus");
                                        if (!e) {
                                            i.style.display = "none";
                                            continue;
                                        }
                                        const c = { columns: [], functions: [] };
                                        for (const e of i.querySelectorAll("td.search-cell")) {
                                            if (!t.test(e.innerText)) continue;
                                            const n = e.getAttribute("data-type");
                                            n && (e.classList.add(`${o}-matching-cell`), c.columns.push(n));
                                        }
                                        let u = n[e];
                                        void 0 === u && (u = { functions: [] });
                                        for (const e of u.functions) {
                                            let n = "";
                                            (t.test(e.name) && (n = `${e.tool}: ${e.product}`),
                                                t.test(e.function) && (n = `${e.tool}: ${e.function}`),
                                                t.test(e.description) && (n = `${e.tool}: ${e.description}`),
                                                n && -1 === c.functions.indexOf(n) && c.functions.push(n));
                                        }
                                        if (0 === c.columns.length && 0 === c.functions.length) {
                                            i.style.display = "none";
                                            continue;
                                        }
                                        ((l = !0), a.push(e));
                                        const h = [];
                                        (c.columns.length > 0 && h.push("Indicated columns<br>"),
                                            c.functions.length > 0 &&
                                                h.push(
                                                    `Gene functions:<ul><li>${c.functions.join("</li><li>")}</li></ul>`,
                                                ),
                                            s(i, ++r, h.join("")));
                                    }
                                    return (
                                        a.length > 0
                                            ? ((0, i.clearSelectedOrfs)(), (0, i.selectOrfsByLoci)(a))
                                            : (0, i.clearSelectedOrfs)(),
                                        a.length > 0
                                    );
                                })(new RegExp(n, "i"), h[0], e.orfs)
                                    ? (r.hide(), h.hide(), g.show())
                                    : (r.show(), h.show(), g.hide()),
                                    f.prop("checked") && (0, i.zoom_to_selection)());
                        })
                        .off("keyup")
                        .on("keyup", (t) => t.stopPropagation()),
                        b[0] instanceof HTMLInputElement && b[0].value.length > 0 && b.trigger("input"),
                        $(`${l} tbody tr`)
                            .off("click")
                            .on("click", function (t) {
                                $(t.target).closest("td").hasClass(`${o}-cell-interactable`) ||
                                    ((v = !1),
                                    i.locusSelectionHandler.call(this, t),
                                    f.prop("checked") && (0, i.zoom_to_selection)(),
                                    setTimeout(() => {
                                        v = !0;
                                    }, 50));
                            }));
                    const x = $(l).find(".seq-copy");
                    (x.each(function () {
                        const e = this.getAttribute("data-locus");
                        if (void 0 === e) return;
                        const n = t.orfs.filter((t) => t.locus_tag === e)[0];
                        this.classList.contains("seq-copy-aa")
                            ? this.setAttribute("data-seq", n.translation)
                            : this.classList.contains("seq-copy-nt") && this.setAttribute("data-seq", n.dna);
                    }),
                        x.on("click", r.copyToClipboard),
                        (0, a.replaceWildcards)(h[0], t));
                };
            },
            82: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.redrawGenericDomains = e.drawGenericDomains = void 0));
                const r = n(86),
                    i = n(677),
                    a = n(220),
                    o = n(3);
                let s;
                const l = { text_height: 14, unique_id: 0, version: "0.0.1" };
                let c = null,
                    u = 25;
                const h = "generic-domain";
                function d(t, e, n, r, i, a, o, s, c, u = !1) {
                    const d = (a + i) * n + 4,
                        g = t.append("g").attr("class", "domain-group");
                    (g
                        .append("text")
                        .text(e.id)
                        .attr("x", 0)
                        .attr("y", d + 0.7 * a)
                        .attr("class", `${h}-orflabel`)
                        .attr("data-locus", e.id),
                        g
                            .append("line")
                            .attr("y1", d + a / 2)
                            .attr("y2", d + a / 2)
                            .attr("x1", s(0))
                            .attr("x2", s(e.seqLength))
                            .attr("class", `${h}-line`),
                        g
                            .append("line")
                            .attr("x1", s(0))
                            .attr("x2", s(0))
                            .attr("y1", d + a / 4)
                            .attr("y2", d + (3 * a) / 4)
                            .attr("class", `${h}-line`),
                        g
                            .append("line")
                            .attr("x1", s(e.seqLength))
                            .attr("x2", s(e.seqLength))
                            .attr("y1", d + a / 4)
                            .attr("y2", d + (3 * a) / 4)
                            .attr("class", `${h}-line`),
                        g
                            .selectAll(`rect.${h}-domain`)
                            .data(e.domains)
                            .enter()
                            .append("rect")
                            .attr("x", (t) => s(t.start))
                            .attr("y", d)
                            .attr("rx", 17)
                            .attr("ry", 17)
                            .attr("width", (t) => s(t.end) - s(t.start))
                            .attr("height", a)
                            .attr("data-id", (t, e) => `${c.name}-details-orf-${r}-${e}-tooltip`)
                            .attr("class", (t) => `${h}-domain ${t.html_class ? t.html_class : "generic-type-other"}`)
                            .attr("stroke-width", (t) => (t.go_terms && t.go_terms.length > 0 ? 3 : 1)),
                        g
                            .selectAll(`text.${h}-text`)
                            .data(e.domains)
                            .enter()
                            .append("text")
                            .text((t) => {
                                if (u) return t.name;
                                const e = s(t.end) - s(t.start),
                                    n = (t.name.match(/[a-z]/g) || []).length,
                                    r = (t.name.match(/[A-Z]/g) || []).length;
                                return e > 16 * r + 18 * (t.name.length - r - n) + 8 * n ? t.name : "...";
                            })
                            .attr("x", (t) => s((t.start + t.end) / 2))
                            .attr("text-anchor", "middle")
                            .attr("y", d + 0.7 * a)
                            .attr("class", (t) => `${h}-text ${t.html_class}`)
                            .attr("font-size", l.text_height)
                            .attr("font-weight", "bold"));
                }
                function g(t, e, n) {
                    ((c = e), (u = n));
                    for (const r of e)
                        $(`.${t}-${r.name}-details`)
                            .off(".firstClick")
                            .one("click.firstClick", () => {
                                f(`${t}-${r.name}-details-svg`, r, n);
                            });
                }
                function f(t, e, n) {
                    if ($(`#${t}`).find(`svg.${h}-svg`).length > 0) return;
                    const s = $("input.domains-expand-full").prop("checked"),
                        c = (0, i.select)(`#${t}`),
                        u = n;
                    let g = $(`#${t}`).width() || 1200;
                    const f = (u + 10) * e.data.length,
                        y = c.append("svg").attr("height", f).attr("width", g).attr("class", `${h}-svg`);
                    let b = 0,
                        x = "";
                    for (const t of e.data) ((b = Math.max(b, t.seqLength || 0)), x.length < t.id.length && (x = t.id));
                    const w = y.append("g"),
                        A = w.append("text").text(x).attr("x", 0).attr("y", 0).attr("id", "dummy-label"),
                        _ = A.node().getComputedTextLength() || 10 * x.length;
                    if (s) {
                        let t = 0;
                        A.attr("class", `${h}-text`);
                        for (const n of e.data)
                            for (const e of n.domains) {
                                A.text(e.name);
                                const n = e.end - e.start,
                                    r = A.node().getComputedTextLength() / n;
                                r > t && (t = r);
                            }
                        const n = t * b;
                        (n > g && (g = n), y.attr("width", g).attr("viewbox", `-1 -1 ${g} ${f}`));
                    }
                    w.remove();
                    const C = (0, r.scaleLinear)()
                            .domain([1, 1.02 * b])
                            .range([_ + 10, g]),
                        S = c.append("div").attr("class", `${h}-svg-singles`),
                        M = u + 8;
                    for (let t = 0; t < e.data.length; t++) {
                        const n = e.data[t],
                            r = l.unique_id++,
                            i = (0, o.locusToFullId)(n.id),
                            f = `${i}-generic-domains`;
                        (d(
                            S.append("svg")
                                .attr("height", M)
                                .attr("width", s ? g : "100%")
                                .attr("viewbox", `-1 -1 ${g} ${M}`)
                                .attr("class", `${h}-svg-single ${f}`),
                            n,
                            0,
                            r,
                            10,
                            u,
                            0,
                            C,
                            e,
                            s,
                        ),
                            d(y, n, t, r, 10, u, 0, C, e, s),
                            $(`#${i}-svgeneorf`)[0].classList.contains("svgene-selected-orf") || $(`.${f}`).hide(),
                            c
                                .append("div")
                                .attr("id", `${e.name}-details-orf-${r}`)
                                .selectAll(`div.${h}-tooltip`)
                                .data(n.domains)
                                .enter()
                                .append("div")
                                .attr("class", `${h}-tooltip`)
                                .attr("id", (t, n) => `${e.name}-details-orf-${r}-${n}-tooltip`)
                                .html((t) => p(t, e.url)),
                            $(`.${h}-tooltip .clipboard-copy`).off("click").click(a.copyToClipboard));
                    }
                    ($(`.${h}-orflabel`)
                        .off("click")
                        .click(function (t) {
                            ($(`#${(0, o.locusToFullId)($(this).attr("data-locus") || "none")}-svgeneorf`).trigger(t),
                                (0, o.zoom_to_selection)());
                        }),
                        (0, i.selectAll)("g.domain-group").data(e.data),
                        $(`.${h}-domain`).click(m),
                        v());
                }
                function p(t, e) {
                    let n = `${t.name}<br>`;
                    if (e.length > 0) {
                        const r = `<a class='external-link' target='_blank' href="${e.replace(
                            "$ACCESSION",
                            t.accession,
                        )}">${t.accession}</a>`;
                        n = t.accession !== t.name ? `${r} - ${t.name}<br>` : `${r}<br>`;
                    }
                    if (
                        ((n += `Location: ${t.start}-${t.end} AA<br>`),
                        (n += `Score: ${t.score}, E-Value: ${t.evalue}<br>`),
                        (n += `${t.description}<br><br>`),
                        t.go_terms && t.go_terms.length > 0)
                    ) {
                        n += "Gene Ontologies:<br>";
                        for (const e of t.go_terms) n += `&nbsp;${e}<br>`;
                        n += "<br>";
                    }
                    return n;
                }
                function m(t) {
                    const e = $(this).attr("data-id");
                    if (void 0 === e) return;
                    const n = $(`#${e}`);
                    if ((s && (clearTimeout(n.data("timeout")), s.hide()), (s = n), "none" !== n.css("display")))
                        return void n.hide();
                    let r = $(this).offset();
                    void 0 === r && (r = { left: 0, top: 0 });
                    let i = setTimeout(() => n.slideUp("fast"), 5e3);
                    n.css("top", r.top + 30)
                        .css("left", r.left + 10)
                        .show()
                        .click(function () {
                            $(this).hide();
                        })
                        .data("timeout", i)
                        .mouseover(() => clearTimeout(n.data("timeout")))
                        .mouseout(() => {
                            (clearTimeout(n.data("timeout")),
                                (i = setTimeout(() => n.slideUp("fast"), 5e3)),
                                n.data("timeout", i));
                        });
                }
                function v(t) {
                    if (t && t.reset) {
                        g(t.anchor, c, u);
                        for (const e of c) {
                            const n = `${t.anchor}-${e.name}-details`;
                            ($(`#${n}-svg`).empty(),
                                $(`.${n}`).hasClass("body-details-header-active") &&
                                    ($(`.${n}`).off(".firstClick"), f(`${n}-svg`, e, u)));
                        }
                    }
                    $("input.domains-selected-only").prop("checked")
                        ? ($(`.${h}-svg`).hide(), $(`.${h}-svg-singles`).show())
                        : ($(`.${h}-svg`).show(), $(`.${h}-svg-singles`).hide());
                }
                ((e.drawGenericDomains = g), (e.redrawGenericDomains = v));
            },
            878: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.createModuleHandlers = e.redrawDomains = e.drawDomains = void 0));
                const r = n(777),
                    i = n(86),
                    a = n(677),
                    o = n(220),
                    s = n(3),
                    l = n(130);
                let c;
                const u = { text_height: 14, unique_id: 0, version: "0.0.1" };
                function h(t, e, n, i, a) {
                    const o = (0, r.path)(),
                        s = e(t.start),
                        l = e(t.end),
                        c = n + i;
                    return (
                        "head" === t.multi_cds
                            ? (o.moveTo(l + a, n),
                              o.lineTo(l, n + i / 3),
                              o.lineTo(l + a, n + (2 * i) / 3),
                              o.lineTo(l, c),
                              o.lineTo(s + a, c),
                              o.arcTo(s, c, s, n, a),
                              o.lineTo(s, n + a),
                              o.arcTo(s, n, l, n, a))
                            : (o.moveTo(l - a, n),
                              o.arcTo(l, n, l, c, a),
                              o.lineTo(l, c - a),
                              o.arcTo(l, c, s, c, a),
                              o.lineTo(s - a, c),
                              o.lineTo(s, n + (2 * i) / 3),
                              o.lineTo(s - a, n + i / 3),
                              o.lineTo(s, n)),
                        o.closePath(),
                        o.toString()
                    );
                }
                function d(t, e, n, i, a, o, s, l) {
                    const c = (o + a) * n + 4,
                        d = t.append("g").attr("class", "domain-group");
                    (d
                        .append("text")
                        .text(e.id)
                        .attr("x", 0)
                        .attr("y", c + 0.7 * o)
                        .attr("class", "jsdomain-orflabel")
                        .attr("data-locus", e.id),
                        d
                            .selectAll("rect.jsdomain-module")
                            .data(e.modules.filter((t) => !t.multi_cds))
                            .enter()
                            .append("rect")
                            .attr("x", (t) => l(t.start))
                            .attr("y", c - 3)
                            .attr("width", (t) => l(t.end) - l(t.start))
                            .attr("height", o + 6)
                            .attr("rx", 5)
                            .attr("class", (t) =>
                                t.complete ? "jsdomain-module" : "jsdomain-module jsdomain-incomplete-module",
                            ),
                        d
                            .selectAll("text.jsdomain-match")
                            .data(e.modules.filter((t) => t.multi_cds))
                            .enter()
                            .append("text")
                            .text((t) => t.match_id)
                            .attr(
                                "x",
                                (t) => l("head" === t.multi_cds ? t.end : t.start) + ("head" === t.multi_cds ? 7 : -5),
                            )
                            .attr("y", (t) => c + (("head" === t.multi_cds ? 2 : 4) * o) / 4 - 2)
                            .attr("text-anchor", "middle")
                            .attr("class", "jsdomain-match"),
                        d
                            .selectAll("path.jsdomain-module")
                            .data(e.modules.filter((t) => t.multi_cds))
                            .enter()
                            .append("path")
                            .attr("d", (t) => h(t, l, c - 3, o + 6, 5))
                            .attr("class", "jsdomain-module"),
                        d
                            .append("line")
                            .attr("y1", c + o / 2)
                            .attr("y2", c + o / 2)
                            .attr("x1", l(0))
                            .attr("x2", l(e.sequence.length))
                            .attr("class", "jsdomain-line"),
                        d
                            .append("line")
                            .attr("x1", l(0))
                            .attr("x2", l(0))
                            .attr("y1", c + o / 4)
                            .attr("y2", c + (3 * o) / 4)
                            .attr("class", "jsdomain-line"),
                        d
                            .append("line")
                            .attr("x1", l(e.sequence.length))
                            .attr("x2", l(e.sequence.length))
                            .attr("y1", c + o / 4)
                            .attr("y2", c + (3 * o) / 4)
                            .attr("class", "jsdomain-line"),
                        d
                            .selectAll("rect.jsdomain-domain")
                            .data(e.domains)
                            .enter()
                            .append("rect")
                            .attr("x", (t) => l(t.start))
                            .attr("y", c)
                            .attr("rx", 17)
                            .attr("ry", 17)
                            .attr("width", (t) => l(t.end) - l(t.start))
                            .attr("height", o)
                            .attr("data-id", (t, e) => `details-orf-${i}-${e}-tooltip`)
                            .attr("class", (t) => `jsdomain-domain ${t.html_class}`)
                            .attr("stroke-width", 1),
                        d
                            .selectAll("text.jsdomain-text")
                            .data(e.domains)
                            .enter()
                            .append("text")
                            .text((t) => t.abbreviation)
                            .attr("x", (t) => l((t.start + t.end) / 2))
                            .attr("text-anchor", "middle")
                            .attr("y", c + 0.7 * o)
                            .attr("class", "jsdomain-text")
                            .attr("font-size", u.text_height)
                            .attr("font-weight", "bold"));
                    const g = d
                        .selectAll("g.jsdomain-module-lid")
                        .data(e.modules.filter((t) => t.complete && !t.multi_cds))
                        .enter()
                        .append("g")
                        .attr("class", "jsdomain-module-lid");
                    (g
                        .append("rect")
                        .attr("x", (t) => l(t.start))
                        .attr("y", c - 3)
                        .attr("width", (t) => l(t.end) - l(t.start))
                        .attr("height", o + 6)
                        .attr("rx", 5)
                        .attr("class", "jsdomain-module-lid-body"),
                        g
                            .append("text")
                            .text((t) => t.monomer || "no prediction")
                            .attr("x", (t) => l((t.start + t.end) / 2))
                            .attr("y", (t) => c + o * (t.iterative ? 0.5 : 0.7))
                            .attr("text-anchor", "middle")
                            .attr("class", "jsdomain-text"));
                    const f = d
                        .selectAll("g.jsdomain-module-lid-fragment")
                        .data(e.modules.filter((t) => t.multi_cds))
                        .enter()
                        .append("g")
                        .attr("class", "jsdomain-module-lid jsdomain-module-lid-fragment");
                    (f
                        .append("path")
                        .attr("d", (t) => h(t, l, c - 3, o + 6, 5))
                        .attr("class", "jsdomain-module-lid-body"),
                        f
                            .append("text")
                            .text((t) => t.monomer || "no prediction")
                            .attr("x", (t) => l((t.start + t.end) / 2))
                            .attr("y", (t) => c + o * (t.iterative ? 0.5 : 0.7))
                            .attr("text-anchor", "middle")
                            .attr("class", "jsdomain-text"));
                    const p = (0, r.path)();
                    (p.moveTo(50, 50),
                        p.arc(0, 50, 50, 0, (Math.PI / 180) * 45, !0),
                        p.moveTo(50, 50),
                        p.lineTo(10, 40));
                    const m = p.toString();
                    g.append("path")
                        .attr("d", m)
                        .attr("stroke", (t) => (t.iterative ? "black" : "none"))
                        .attr("stroke-width", "5px")
                        .attr("fill", "none")
                        .attr("stroke-linecap", "round")
                        .attr(
                            "transform",
                            (t) => `translate(${l((t.start + t.end) / 2)}, ${c + 0.6 * o}) scale(0.1, 0.1)`,
                        );
                }
                function g(t, e) {
                    let n = `${t.type}<br>Location: ${t.start}-${t.end} AA<br>`;
                    if (
                        (t.napdoslink.length > 0 &&
                            (n += `<a href="${t.napdoslink}" target="_blank">Analyze with NaPDoS</a><br>`),
                        (n += `<a href="${t.blastlink}" target="_blank">Run BlastP on this domain</a><br>`),
                        t.predictions.length > 0)
                    ) {
                        n += "<dl><dt>Predictions:</dt>";
                        for (const e of t.predictions) n += `<dd>-${e[0]}: ${e[1]}</dd>`;
                        n += "</dl>";
                    }
                    return (
                        (n += `AA sequence: ${(0, o.clipboardCopyConstruct)(t.sequence)}<br>`),
                        t.dna_sequence &&
                            (n += `Nucleotide sequence: ${(0, o.clipboardCopyConstruct)(t.dna_sequence)}<br>`),
                        (n = (0, l.replaceWildcardsInText)(n, t)),
                        n
                    );
                }
                function f(t) {
                    const e = $(this).attr("data-id");
                    if (void 0 === e) return;
                    const n = $(`#${e}`);
                    if ((c && (clearTimeout(n.data("timeout")), c.hide()), (c = n), "none" !== n.css("display")))
                        return void n.hide();
                    let r = $(this).offset();
                    void 0 === r && (r = { left: 0, top: 0 });
                    let i = setTimeout(() => n.slideUp("fast"), 5e3);
                    n.css("top", r.top + 30)
                        .css("left", r.left + 10)
                        .show()
                        .click(function () {
                            $(this).hide();
                        })
                        .data("timeout", i)
                        .mouseover(() => clearTimeout(n.data("timeout")))
                        .mouseout(() => {
                            (clearTimeout(n.data("timeout")),
                                (i = setTimeout(() => n.slideUp("fast"), 5e3)),
                                n.data("timeout", i));
                        });
                }
                function p(t) {
                    if (($(this).hide(), t.clientX && t.clientY)) {
                        const e = document.elementFromPoint(t.clientX, t.clientY);
                        e && $(e).click();
                    }
                    $(this).show();
                }
                function m() {
                    ($("input.domains-selected-only").prop("checked")
                        ? ($(".jsdomain-svg").hide(), $(".jsdomain-svg-singles").show())
                        : ($(".jsdomain-svg").show(), $(".jsdomain-svg-singles").hide()),
                        $("input.show-module-domains").prop("checked")
                            ? $(".jsdomain-module-lid").hide()
                            : $(".jsdomain-module-lid").show());
                }
                ((e.drawDomains = function (t, e, n) {
                    $(`.body-details-header.${t}-nrps_pks`)
                        .off(".firstClick")
                        .one("click.firstClick", () => {
                            !(function (t, e, n) {
                                const r = (0, a.select)(`#${t}`),
                                    l = n,
                                    c = $(`#${t}`).width() || 700;
                                (r.selectAll("svg.jsdomain-svg").remove(),
                                    r.selectAll("svg.jsdomain-svg-single").remove());
                                const h = (l + 10) * e.orfs.length,
                                    v = r
                                        .append("svg")
                                        .attr("height", h)
                                        .attr("width", "100%")
                                        .attr("viewbox", `-1 -1 ${c} ${h}`)
                                        .attr("class", "jsdomain-svg");
                                let y = 0,
                                    b = "";
                                for (const t of e.orfs)
                                    ((y = Math.max(y, t.sequence.length)), b.length < t.id.length && (b = t.id));
                                const x = v
                                        .append("g")
                                        .append("text")
                                        .text(b)
                                        .attr("x", 0)
                                        .attr("y", 0)
                                        .attr("id", "dummy-label"),
                                    w = x.node().getComputedTextLength() || 15 * b.length;
                                x.remove();
                                const A = (0, i.scaleLinear)()
                                        .domain([1, 1.02 * y])
                                        .range([w + 10, c]),
                                    _ = r.append("div").attr("class", "jsdomain-svg-singles"),
                                    C = l + 8;
                                for (let t = 0; t < e.orfs.length; t++) {
                                    const n = e.orfs[t],
                                        i = u.unique_id++;
                                    (d(
                                        _.append("svg")
                                            .attr("height", C)
                                            .attr("width", "100%")
                                            .attr("viewbox", `-1 -1 ${c} ${C}`)
                                            .attr("class", "jsdomain-svg-single")
                                            .attr("id", `${(0, s.locusToFullId)(n.id)}-domains`),
                                        n,
                                        0,
                                        i,
                                        10,
                                        l,
                                        0,
                                        A,
                                    ),
                                        d(v, n, t, i, 10, l, 0, A),
                                        r
                                            .append("div")
                                            .attr("id", `details-orf-${i}`)
                                            .selectAll("div.jsdomain-tooltip")
                                            .data(n.domains)
                                            .enter()
                                            .append("div")
                                            .attr("class", "jsdomain-tooltip")
                                            .attr("id", (t, e) => `details-orf-${i}-${e}-tooltip`)
                                            .html((t) => g(t)),
                                        $(".jsdomain-tooltip .clipboard-copy").off("click").click(o.copyToClipboard));
                                }
                                ($(".jsdomain-orflabel")
                                    .off("click")
                                    .click(function (t) {
                                        $(
                                            `#${(0, s.locusToFullId)($(this).attr("data-locus") || "none")}-svgeneorf`,
                                        ).trigger(t);
                                    }),
                                    (0, a.selectAll)("g.domain-group").data(e.orfs),
                                    $(".jsdomain-domain").click(f),
                                    $(".jsdomain-module-lid").click(p),
                                    $(".jsdomain-svg-single").hide(),
                                    $(".svgene-selected-orf").each(function () {
                                        const t = ($(this).attr("id") || "invalid-id").replace(
                                            "-svgeneorf",
                                            "-domains",
                                        );
                                        $(`#${t}`).show();
                                    }),
                                    m());
                            })(`${t}-details-svg`, e, n);
                        });
                }),
                    (e.redrawDomains = m),
                    (e.createModuleHandlers = function () {
                        $("input.show-module-domains").change(function () {
                            ($("input.show-module-domains").prop("checked", $(this).prop("checked")), m());
                        });
                    }));
            },
            681: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.createRecordOverviews = void 0));
                const r = n(86),
                    i = n(677);
                function a(t, e, n) {
                    if (0 === e.length) return;
                    const a = (0, i.select)(`#record-minimap-${t}`);
                    a.attr("width", 810).attr("height", 40);
                    const o = (0, r.scaleLinear)().domain([0, n]).range([0, 800]);
                    (a
                        .append("line")
                        .attr("x1", 0)
                        .attr("y1", 20)
                        .attr("x2", 800)
                        .attr("y2", 20)
                        .attr("class", "minimap-record-end centerline"),
                        a
                            .append("line")
                            .attr("x1", 0)
                            .attr("y1", 0)
                            .attr("x2", 0)
                            .attr("y2", 40)
                            .attr("class", "minimap-record-end centerline"),
                        a
                            .append("line")
                            .attr("x1", 800)
                            .attr("y1", 0)
                            .attr("x2", 800)
                            .attr("y2", 40)
                            .attr("class", "centerline"),
                        a
                            .selectAll("line.minimap-region-label-line")
                            .data(e)
                            .enter()
                            .append("line")
                            .attr("class", "minimap-region-label-line centerline")
                            .attr("x1", (t) => o((t.start + t.end) / 2))
                            .attr("y1", (t) => (parseInt(t.anchor.split("c")[1], 10) % 2 == 0 ? 25 : 15))
                            .attr("x2", (t) => o((t.start + t.end) / 2))
                            .attr("y2", (t) => (parseInt(t.anchor.split("c")[1], 10) % 2 == 0 ? 30 : 10)),
                        a
                            .selectAll("rect.minimap-region")
                            .data(e)
                            .enter()
                            .append("rect")
                            .attr("class", (t) => `minimap-region ${t.cssClass}`)
                            .attr("y", 15)
                            .attr("x", (t) => o(t.start))
                            .attr("width", (t) => o(t.end - t.start))
                            .attr("height", 10)
                            .style("stroke", "black"),
                        a
                            .selectAll("text.minimap-region-label")
                            .data(e)
                            .enter()
                            .append("text")
                            .attr("x", (t) => o((t.start + t.end) / 2))
                            .attr("y", 20)
                            .attr("font-size", "xx-small")
                            .attr("class", "minimap-region-label")
                            .attr("text-anchor", "middle")
                            .attr("dy", (t) => (parseInt(t.anchor.split("c")[1], 10) % 2 == 0 ? "2" : "-1.2") + "em")
                            .text((t) => t.anchor.split("c")[1]));
                }
                e.createRecordOverviews = function (t) {
                    let e = 0;
                    for (const n of t) ((e += 1), 0 !== n.regions.length && a(e, n.regions, n.length));
                };
            },
            625: function (t, e, n) {
                "use strict";
                var r =
                    (this && this.__importDefault) ||
                    function (t) {
                        return t && t.__esModule ? t : { default: t };
                    };
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.drawStructures = e.drawStructure = void 0));
                const i = r(n(261));
                function a() {
                    const t = this,
                        e = { height: 200, padding: 30, width: ($(t).parent().width() || 270) - 20 },
                        n = $(t).parent(),
                        r = $(t).is(":visible");
                    r || (n.detach(t.id), $("body").append(t));
                    const a = new i.default.Drawer(e);
                    (i.default.parse(t.getAttribute("data-smiles"), (e) => a.draw(e, t, "light", !1)),
                        r || n.prepend(t),
                        $(this).click(function (t) {
                            $(`#${$(this).attr("id")}-modal`).show();
                            const e = $(`#${$(this).attr("id")}-modal`)
                                .children(".modal-content")
                                .children(".smiles-canvas-modal")
                                .first();
                            if (void 0 === e) return;
                            const n = { height: e.height() || 500, padding: 30, width: (e.width() || 700) + 80 },
                                r = e.parent().parent().hasClass("expanded");
                            r || e.parent().css("display", "block");
                            const a = new i.default.Drawer(n),
                                o = e[0];
                            (i.default.parse(o.getAttribute("data-smiles"), (t) => a.draw(t, o, "light", !1)),
                                r || $(e).parent().removeAttr("style"),
                                $(`#${$(this).attr("id")}-modal`)
                                    .children(".modal-content")
                                    .css("top", $(this).position().top - n.height / 2 + "px"));
                        }));
                }
                ((e.drawStructure = a),
                    (e.drawStructures = function () {
                        ($(".smiles-canvas").each(a),
                            $(".modal-container")
                                .off("click")
                                .click(() => $(".modal-container").hide()));
                    }));
            },
            24: (t, e, n) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }), (e.drawBindingSites = void 0));
                const r = n(677),
                    i = n(3),
                    a = 4.5,
                    o = 2.25,
                    s = 40;
                function l(t, e, n, r, i) {
                    const a = i / 2;
                    (t
                        .append("line")
                        .attr("x1", e)
                        .attr("x2", n)
                        .attr("y1", r)
                        .attr("y2", r)
                        .attr("class", "tfbs-arrow-base"),
                        e > n && (i *= -1),
                        t
                            .append("path")
                            .attr("d", `M ${n} ${r} L ${n - i} ${r + a} L ${n - i} ${r - a} z`)
                            .attr("class", "tfbs-arrow-head"));
                }
                e.drawBindingSites = function (t, e) {
                    if (!e) return;
                    let n = 0,
                        c = 0,
                        u = 0,
                        h = 0;
                    for (const t of e)
                        "weak" !== t.confidence
                            ? (t.left && t.left.name.length > n && (n = t.left.name.length),
                              t.right && t.right.name.length > c && (c = t.right.name.length))
                            : (t.left && t.left.name.length > u && (u = t.left.name.length),
                              t.right && t.right.name.length > h && (h = t.right.name.length));
                    const d = Math.max(100, 12 * Math.max(n, c)),
                        g = d + s + 9 * (e[0].presequence.length + 2),
                        f = g + 300;
                    for (const n of e) {
                        const e = `tfbs-${t}-${n.start}-${n.name}`,
                            i = $(`#${e}`);
                        if ($(i).children().length > 0) continue;
                        const c = 300 / n.sequence.length,
                            u = 9 * n.presequence.length + a,
                            h = 9 * n.postsequence.length + a,
                            p = g - u - a,
                            m = f + h + a,
                            v = (0, r.select)(`#${e}`)
                                .append("svg")
                                .attr("width", m + s + d)
                                .attr("height", 99)
                                .attr("id", e)
                                .attr("class", "tfbs-svg");
                        v.append("text")
                            .attr("x", p)
                            .attr("y", 40)
                            .attr("textLength", u)
                            .attr("class", "tfbs-sequence tfbs-hit-text-nearby")
                            .text(n.presequence);
                        const y = v.append("g");
                        v.append("text")
                            .attr("x", f + a)
                            .attr("y", 40)
                            .attr("textLength", u)
                            .attr("class", "tfbs-sequence tfbs-hit-text-nearby")
                            .text(n.postsequence);
                        const b = v.append("g"),
                            x = v.append("g");
                        for (let t = 0; t < n.target.length; t++) {
                            const e = g + t * c + c / 4;
                            (y
                                .append("text")
                                .attr("x", e)
                                .attr("y", 40)
                                .attr("class", "tfbs-sequence")
                                .text(n.sequence[t]),
                                b
                                    .append("text")
                                    .attr("x", e)
                                    .attr("y", 58)
                                    .attr("class", "tfbs-sequence tfbs-hit-text-nearby")
                                    .text(n.matches[t] ? "|" : " "),
                                x
                                    .append("text")
                                    .attr("x", e)
                                    .attr("y", 74)
                                    .attr(
                                        "class",
                                        "tfbs-sequence " + ("N" === n.target[t] ? "tfbs-hit-text-nearby" : ""),
                                    )
                                    .text(n.target[t]));
                        }
                        const w = v.append("g").attr("class", "tfbs-boundary");
                        (w
                            .append("text")
                            .text(n.start.toLocaleString())
                            .attr("text-anchor", "middle")
                            .attr("x", g - o)
                            .attr("y", 20)
                            .attr("class", "tfbs-coordinate"),
                            w
                                .append("line")
                                .attr("x1", g - o)
                                .attr("x2", g - o)
                                .attr("y1", 25)
                                .attr("y2", 79)
                                .attr("class", "tfbs-line"),
                            w
                                .append("text")
                                .text(n.end.toLocaleString())
                                .attr("text-anchor", "middle")
                                .attr("x", f + o)
                                .attr("y", 20)
                                .attr("class", "tfbs-coordinate"),
                            w
                                .append("line")
                                .attr("x1", f + o)
                                .attr("x2", f + o)
                                .attr("y1", 25)
                                .attr("y2", 79)
                                .attr("class", "tfbs-line"));
                        const A = (t, e, n = 0) => {
                                const r = v.append("g").attr("class", "tfbs-nearby");
                                for (let n = 0; n < 3; n++)
                                    r.append("circle")
                                        .attr("cx", t + (5 + 10 * n) * e)
                                        .attr("cy", 34)
                                        .attr("r", 2)
                                        .attr("class", "tfbs-ellipsis");
                                if (!(n <= 0)) {
                                    for (let n = -1; n < 2; n += 2) l(r, t + 15 * e, t + 15 * e + 30 * n, 74, 8);
                                    r.append("text")
                                        .text(n.toLocaleString())
                                        .attr("x", t + 15 * e)
                                        .attr("y", 66)
                                        .attr("text-anchor", "middle")
                                        .attr("class", "tfbs-hit-text-nearby");
                                }
                            },
                            _ = (t, e, n, r, i, a, o) => {
                                let s = 0.05 * d;
                                i < 0 && (s *= -1);
                                let l = "";
                                l = r
                                    ? `M ${n},45\n                        L ${n + i},45\n                        L ${
                                          n + i
                                      },25\n                        L ${n},25`
                                    : `M ${n},45\n                        L ${
                                          n + i - s
                                      },45\n                        L ${n + i},35\n                        L ${
                                          n + i - s
                                      },25\n                        L ${n},25`;
                                const c = t.append("g").attr("class", "tfbs-gene");
                                (c
                                    .append("text")
                                    .attr("x", o)
                                    .attr("y", 40)
                                    .text(e.name)
                                    .attr("text-anchor", a)
                                    .attr("class", "serif tfbs-gene-label")
                                    .attr("data-locus", e.name),
                                    c.append("path").attr("class", "tfbs-line").attr("d", l).attr("fill", "none"));
                            };
                        if (n.contained_by_left)
                            (_(v, n.left, p - s, 1 === n.left.strand, -d, "start", p - s - 0.9 * d),
                                A(p, -1, n.start - n.left.location),
                                _(v, n.right, m + s, -1 === n.left.strand, d, "end", m + s + 0.9 * d),
                                A(m, 1, n.right.location - n.start + n.sequence.length));
                        else {
                            if (n.left) {
                                const t = n.start - n.left.location,
                                    e = p - d - s;
                                let r = d;
                                (t <= 0 && (r = g - (t + 0.2) * c - e),
                                    _(v, n.left, e, -1 === n.left.strand, r, "end", p - s - 0.1 * d),
                                    A(p, -1, t));
                            }
                            if (n.mid) {
                                const t = g + (n.mid.location - n.start) * c,
                                    e = t + n.mid.length * c,
                                    r = 0.05 * (e - t);
                                let i = "";
                                ((i =
                                    -1 !== n.mid.strand
                                        ? `M ${t},45\n                            L ${
                                              e - r
                                          },45\n                            L ${e},35\n                            L ${
                                              e - r
                                          },25\n                            L ${t},25\n                            Z`
                                        : `M ${e},45\n                            L ${
                                              t + r
                                          },45\n                            L ${t},35\n                            L ${
                                              t + r
                                          },25\n                            L ${e},25\n                            Z`),
                                    v
                                        .append("g")
                                        .attr("class", "tfbs-gene")
                                        .append("path")
                                        .attr("class", "tfbs-line")
                                        .attr("d", i)
                                        .attr("fill", "none"));
                            }
                            if (n.right) {
                                const t = m + d + s;
                                let e = -d;
                                (n.right.location < n.end && (e = f + (n.right.location - n.end) * c - t),
                                    _(v, n.right, t, 1 === n.right.strand, e, "start", m + s + 0.1 * d),
                                    A(m, 1, n.right.location - n.end));
                            }
                        }
                        const C = (m - p) / 3,
                            S = p + C,
                            M = m - C,
                            T = v.append("g").attr("class", "tfbs-line tfbs-strand");
                        -1 !== n.strand ? l(T, S, M, 89, 10) : l(T, M, S, 89, 10);
                    }
                    $(".tfbs-gene-label")
                        .off("click")
                        .on("click", function () {
                            const t = this.getAttribute("data-locus");
                            t && (0, i.selectOrfsByLoci)([t]);
                        });
                };
            },
            3: function (t, e, n) {
                "use strict";
                var r =
                    (this && this.__awaiter) ||
                    function (t, e, n, r) {
                        return new (n || (n = Promise))(function (i, a) {
                            function o(t) {
                                try {
                                    l(r.next(t));
                                } catch (t) {
                                    a(t);
                                }
                            }
                            function s(t) {
                                try {
                                    l(r.throw(t));
                                } catch (t) {
                                    a(t);
                                }
                            }
                            function l(t) {
                                var e;
                                t.done
                                    ? i(t.value)
                                    : ((e = t.value),
                                      e instanceof n
                                          ? e
                                          : new n(function (t) {
                                                t(e);
                                            })).then(o, s);
                            }
                            l((r = r.apply(t, e || [])).next());
                        });
                    };
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.resetZoom =
                        e.resetView =
                        e.zoom_to_selection =
                        e.locusSelectionHandler =
                        e.tag_to_id =
                        e.drawRegion =
                        e.clearSelectedOrfs =
                        e.locusToFullId =
                        e.selectOrfsByLoci =
                            void 0));
                const i = n(216),
                    a = n(305),
                    o = n(86),
                    s = n(677);
                n(473);
                const l = n(220),
                    c = n(26),
                    u = n(130);
                let h = 100;
                const d = 14,
                    g = "svgene-selected-orf";
                function f(t, e = !1) {
                    if (!q) return;
                    if (t.length < 1) return;
                    let n = $(`#${p(t[0])}-svgeneorf`);
                    for (const e of t.slice(1)) n = n.add(`#${p(e)}-svgeneorf`);
                    e ? w(n) : C(n, t.length > 1);
                }
                function p(t) {
                    return q ? `u${j - 1}-region${q.idx}-${x(t)}` : `${x(t)}`;
                }
                function m(t) {
                    const e = z + d + F,
                        n = z + d + h - F,
                        r = z + d + h / 2;
                    if (1 === t.strand) {
                        const i = H(t.start),
                            a = Math.max(H(t.end) - 2 * F, i);
                        return `${i},${e} ${a},${e} ${H(t.end)},${r} ${a},${n} ${i},${n}`;
                    }
                    if (-1 === t.strand) {
                        const i = H(t.start),
                            a = H(t.end),
                            o = Math.min(H(t.start) + 2 * F, a);
                        return `${i},${r} ${o},${e} ${a},${e} ${a},${n} ${o},${n}`;
                    }
                    return `${t.start},${e} ${t.end},${e} ${t.end},${n} ${t.start},${n}`;
                }
                function v(t, e, n, r) {
                    const i = n + d + e,
                        a = n + 2 * d + e - r,
                        o = Math.floor(H((t.start + t.end) / 2));
                    return `${o},${i} ${o - 5},${a} ${o + 5},${a}`;
                }
                function y(t, e, n) {
                    return Math.max(t, Math.min(e, n));
                }
                function b(t, e) {
                    return ("biosynthetic" !== t.type && "biosynthetic" !== e.type) ||
                        ("biosynthetic" === t.type && "biosynthetic" === e.type)
                        ? t.start - e.start
                        : "biosynthetic" === t.type
                          ? 1
                          : -1;
                }
                function x(t) {
                    return t.replace(/(:|\.)/g, "-").replace(/-svgeneorf/g, "_orf");
                }
                function w(t) {
                    (t[0].classList.contains(g)
                        ? (L(t), 0 === $(`.${g}`).length && (N(), R()))
                        : (0 === $(`.${g}`).length && L(), N(t)),
                        M(t, "cds"));
                }
                function A(t) {
                    ($(".legend-selected").removeClass("legend-selected"),
                        t.ctrlKey || t.metaKey ? w($(this)) : C($(this)));
                }
                function _(t) {
                    return r(this, void 0, void 0, function* () {
                        try {
                            const t = yield fetch(
                                    "https://antismash-db.secondarymetabolites.org/api/jobs/clusterblast",
                                    {
                                        body: JSON.stringify({
                                            name: $(this).attr("data-locus"),
                                            sequence: $(this).attr("data-seq"),
                                        }),
                                        headers: { "Content-Type": "application/json" },
                                        method: "POST",
                                        mode: "cors",
                                    },
                                ),
                                e = yield t.json();
                            window.open("https://antismash-db.secondarymetabolites.org/job/" + e.id, "_blank");
                        } catch (t) {
                            alert("The antiSMASH database is not currently accessible: " + t);
                        }
                    });
                }
                function C(t, e = !1) {
                    if (!q) return;
                    const n = (0, s.selectAll)(t.toArray()),
                        r = n.datum();
                    if (!e) {
                        const t = $(`.focus-panel-content-${q.anchor}`);
                        (t.html(r.description).find(".collapser").click(c.toggleCollapserHandler),
                            (0, u.replaceWildcards)(t[0], q),
                            $(".clipboard-copy", t).off("click").click(l.copyToClipboard),
                            $(".asdb-linkout", t).off("click").click(_));
                    }
                    (n.classed(g) && 1 === $(`.svgene-orf.${g}`).length ? N() : (L(), N(t)), M(t, "cds"));
                }
                function S(t) {
                    (0, c.toggleCollapser)($(".collapser-level-candidate.expanded"));
                    const e = t.product.split(":")[0].split(" ")[1];
                    (0, c.toggleCollapser)($(`.collapser-level-candidate.collapser-target-CC${e}`));
                }
                function M(t, e) {
                    const n = (0, s.selectAll)(t.toArray()).datum();
                    (0, c.toggleCollapser)($(`.collapser-target-${x(n.locus_tag)}.collapser-level-${e}`));
                }
                function T(t, e) {
                    ($(".legend-selected").removeClass("legend-selected"),
                        $(".svgene-orf").each(function (n) {
                            const r = (0, s.selectAll)($(this).toArray()).datum();
                            t <= r.start &&
                                r.end <= e &&
                                (N($(this)), M($(this), "candidate"), M($(this), "protocluster"));
                        }));
                }
                function k(t, e) {
                    const n = e ? "1" : "0.5",
                        r = (0, s.selectAll)(t.toArray()),
                        i = q && t.length / 2 === q.orfs.length;
                    (i && R(),
                        r
                            .attr("opacity", n)
                            .classed(g, e)
                            .each((t) => {
                                (0, s.selectAll)(`#svgene-minimap-orf-${x(t.locus_tag)}`)
                                    .attr("opacity", n)
                                    .classed(g, e);
                                const r = p(t.locus_tag);
                                (e
                                    ? ($(`#${r}-domains`).show(), $(`.${r}-generic-domains`).show())
                                    : ($(`#${r}-domains`).hide(), $(`.${r}-generic-domains`).hide()),
                                    i ||
                                        (function (t, e, n = !1) {
                                            t &&
                                                (n &&
                                                    (t = (function (t) {
                                                        return t.split(/-/)[2];
                                                    })(t)),
                                                e
                                                    ? $(`.cds-selected-marker-${t}`).addClass("active")
                                                    : $(`.cds-selected-marker-${t}`).removeClass("active"));
                                        })(t.locus_tag, e));
                            }));
                }
                function R() {
                    $(".cds-selected-marker").removeClass("active");
                }
                function N(t) {
                    (void 0 === t && (t = $(".svgene-orf, .svgene-minimap-orf")), k(t, !0));
                }
                function L(t) {
                    (void 0 === t && (t = $(`.${g}`)), k(t, !1));
                }
                function P(t) {
                    const e = $(this).attr("data-id");
                    if (void 0 === e) return;
                    let n = ".contains-tta-codon";
                    if (
                        ((n =
                            "legend-tta-codon" === e
                                ? ".contains-tta-codon"
                                : "legend-resistance" === e
                                  ? ".svgene-resistance-orf"
                                  : `.${e.replace("legend-", "svgene-")}`),
                        $(this).hasClass("legend-selected"))
                    )
                        return (
                            $(this).removeClass("legend-selected"),
                            L($(n)),
                            void (0 === $(".legend-selected").length && N())
                        );
                    (((t.ctrlKey || t.metaKey) && 0 !== $(".legend-selected").length) ||
                        ($(".legend-selected").removeClass("legend-selected"), L()),
                        $(this).addClass("legend-selected"),
                        N($(n)));
                }
                function D(t) {
                    if (!q) return;
                    let e = -1,
                        n = -1;
                    ($(`.${g}`).each(function (t) {
                        const r = (0, s.selectAll)($(this).toArray()).datum();
                        ((e = -1 === e ? Math.min(r.start, r.end) : Math.min(e, r.start, r.end)),
                            (n = -1 === n ? Math.max(r.start, r.end) : Math.max(n, r.start, r.end)));
                    }),
                        -1 !== e && -1 !== n && I(Math.max(q.start, e - 3e3), Math.min(n + 3e3, q.end)));
                }
                function I(t, e, n) {
                    if (!q) return;
                    let r = 0;
                    (void 0 === n && ((n = !1), (r = 500)), (t -= 1), (e += 1), H.domain([t, e]));
                    const i = t + (e - t) / 2;
                    ((0, s.selectAll)(".svgene-orf,.svgene-orf-bg")
                        .transition()
                        .duration(r)
                        .attr("points", (t) => m(t)),
                        (0, s.selectAll)(".svgene-locustag *")
                            .transition()
                            .duration(r)
                            .attr("x", (t) => (t.start < i ? H(t.start) : H(t.end)))
                            .attr("text-anchor", (t) => (t.start < i ? "start" : "end")),
                        (0, s.selectAll)(".cluster-background")
                            .transition()
                            .duration(r)
                            .attr("width", (t) => H(t.neighbouring_end) - H(t.neighbouring_start))
                            .attr("x", (t) => H(t.neighbouring_start)),
                        (0, s.selectAll)(".cluster-core")
                            .transition()
                            .duration(r)
                            .attr("width", (t) => H(t.end) - H(t.start))
                            .attr("x", (t) => H(t.start)),
                        (0, s.selectAll)(".cluster-line")
                            .transition()
                            .duration(r)
                            .attr("x1", (t) => H(t.neighbouring_start))
                            .attr("x2", (t) => H(t.neighbouring_end)),
                        (0, s.selectAll)(".clusterlabel")
                            .transition()
                            .duration(r)
                            .attr("x", (n) => H((Math.max(n.start, t + 1) + Math.min(n.end, e - 1)) / 2)),
                        (0, s.selectAll)(".svgene-resistance")
                            .filter((t) => t && !0 === t.resistance)
                            .transition()
                            .duration(r)
                            .attr("x", (t) => H(t.start))
                            .attr("width", (t) => H(t.end) - H(t.start)),
                        (0, s.selectAll)(".svgene-tta-codon")
                            .filter((t) => void 0 !== t)
                            .transition()
                            .duration(r)
                            .attr("points", (t) => v(t, h, z, F)));
                    let a = (0, s.selectAll)(".svgene-binding-site");
                    if (
                        ((a = a
                            .filter((t) => void 0 !== t)
                            .each(function (t) {
                                const e = H(t.loc + t.len / 2);
                                ((0, s.select)(this)
                                    .selectAll("line")
                                    .transition()
                                    .duration(r)
                                    .attr("x1", e)
                                    .attr("x2", e),
                                    (0, s.select)(this).selectAll("circle").transition().duration(r).attr("cx", e));
                            })),
                        null !== B && (0, s.select)(".svgene-axis").transition().duration(r).call(B),
                        !n)
                    ) {
                        0 === $(`.${g}`).length && N();
                        const n = O(t),
                            i = O(e);
                        ((0, s.select)(`#svgene-minimap-grabber-a-${q.anchor}`).transition().duration(r).attr("x", n),
                            (0, s.select)(`#svgene-minimap-grabber-b-${q.anchor}`)
                                .transition()
                                .duration(r)
                                .attr("x", i),
                            (0, s.select)(`#svgene-minimap-window-${q.anchor}`)
                                .transition()
                                .duration(r)
                                .attr("x", n)
                                .attr("width", i - n));
                    }
                }
                function E() {
                    null !== q &&
                        ($(".legend-selected").removeClass("legend-selected"),
                        N(),
                        I(q.start, q.end),
                        (0, c.toggleCollapser)($(".collapser-level-cds.expanded")));
                }
                ((e.selectOrfsByLoci = f),
                    (e.locusToFullId = p),
                    (e.clearSelectedOrfs = function () {
                        N();
                    }),
                    (e.drawRegion = function (t, e, n, r, l) {
                        (q &&
                            ((0, s.select)(`#${q.anchor}-svg`).selectAll("*").remove(),
                            $(".legend-selected").removeClass("legend-selected")),
                            (q = e));
                        const c = e;
                        if (null === q) return;
                        (void 0 === r && (r = c.start), void 0 === l && (l = c.end));
                        let u = !1;
                        for (const t of c.orfs)
                            if (t.resistance) {
                                u = !0;
                                break;
                            }
                        const f =
                                12 *
                                Math.max.apply(
                                    Math,
                                    c.clusters.map((t) => t.height),
                                ),
                            p = u || c.sites.ttaCodons.length > 0 ? Math.floor((2 * n) / 3) : 0,
                            w = n + 2 * d + f + 40 + 20 + p,
                            _ = (c.idx, (0, s.select)(`#${t}`)),
                            C = $(`#${t}`).parent().width() || 700;
                        (_.selectAll("svg").remove(), _.selectAll("div").remove());
                        const M = _.append("svg")
                                .attr("height", w)
                                .attr("width", "100%")
                                .attr("viewbox", `-1 0 ${C} ${w}`),
                            k = [],
                            R = [],
                            N = c.sites;
                        (k.push.apply(k, c.orfs.sort(b)), R.push.apply(R, c.clusters ? c.clusters : []));
                        const U = j++;
                        ((F = n / 10),
                            (H = (0, o.scaleLinear)().domain([r, l]).range([0, C])),
                            (function (t, e, n, r, i, a, o, l, c, u) {
                                if (null === q) return;
                                const f = q.idx;
                                ((z =
                                    12 *
                                        Math.max.apply(
                                            Math,
                                            n.map((t) => t.height),
                                        ) +
                                    d),
                                    (h = a));
                                const p = z + d + a / 2;
                                t.append("line")
                                    .attr("x1", 0)
                                    .attr("y1", p)
                                    .attr("x2", o)
                                    .attr("y2", p)
                                    .attr("class", "centerline");
                                const y = t
                                    .selectAll("g.cluster-bar-group")
                                    .data(n)
                                    .enter()
                                    .append("g")
                                    .attr("class", (t) =>
                                        "candidatecluster" === t.kind
                                            ? `candidate-${t.product.split(" ")[2].replace("chemical_", "")}`
                                            : "subregion" === t.kind && 0 === t.prefix.length
                                              ? `svgene-border-${t.tool}`
                                              : "",
                                    )
                                    .on("click", (t) => {
                                        (($(`.${g}`).length !== e.length && (s.event.ctrlKey || s.event.metaKey)) ||
                                            L(),
                                            T(t.neighbouring_start, t.neighbouring_end),
                                            S(t));
                                    })
                                    .on("dblclick", (t) => {
                                        (T(t.neighbouring_start, t.neighbouring_end),
                                            S(t),
                                            I(t.neighbouring_start, t.neighbouring_end));
                                    });
                                (y
                                    .append("rect")
                                    .attr("width", (t) => H(t.neighbouring_end) - H(t.neighbouring_start))
                                    .attr("height", 10)
                                    .attr("x", (t) => H(t.neighbouring_start))
                                    .attr("y", (t) => 12 * t.height + l)
                                    .attr("opacity", "0.5")
                                    .attr("class", (t) =>
                                        "subregion" === t.kind
                                            ? "cluster-background"
                                            : `cluster-background ${t.product} ${t.category}`,
                                    )
                                    .style("stroke-width", "0"),
                                    y
                                        .append("line")
                                        .attr("x1", (t) => H(t.neighbouring_start))
                                        .attr("y1", (t) => 12 * t.height + l + 5)
                                        .attr("x2", (t) => H(t.neighbouring_end))
                                        .attr("y2", (t) => 12 * t.height + l + 5)
                                        .attr("class", "cluster-line"),
                                    y
                                        .append("rect")
                                        .attr("width", (t) => H(t.end) - H(t.start))
                                        .attr("height", 10)
                                        .attr("x", (t) => H(t.start))
                                        .attr("y", (t) => 12 * t.height + l)
                                        .attr("class", (t) =>
                                            "subregion" === t.kind
                                                ? t.prefix
                                                    ? "cluster-core"
                                                    : `cluster-core svgene-border-${t.tool}`
                                                : `cluster-core ${t.product} ${t.category}`,
                                        )
                                        .style("stroke", "black"),
                                    y
                                        .append("text")
                                        .attr("x", (t) => H((t.start + t.end) / 2))
                                        .attr("y", (t) =>
                                            "protocluster" === t.kind
                                                ? 12 * (t.height - 1) - 2 + 10 + l
                                                : 12 * t.height - 2 + 10 + l,
                                        )
                                        .style("font-size", "xx-small")
                                        .attr("class", "clusterlabel")
                                        .attr("text-anchor", "middle")
                                        .style("pointer-events", "none")
                                        .text((t) => t.prefix + t.product.replace("_", " ")));
                                const b = p - a / 2 - 4.5,
                                    w = t
                                        .selectAll("g.svgene-binding-site")
                                        .data(r.bindingSites)
                                        .enter()
                                        .append("g")
                                        .attr("class", "svgene-binding-site");
                                (w
                                    .append("line")
                                    .attr("x1", (t) => H(t.loc + t.len / 2))
                                    .attr("x2", (t) => H(t.loc + t.len / 2))
                                    .attr("y1", p)
                                    .attr("y2", b),
                                    w
                                        .append("circle")
                                        .attr("cx", (t) => H(t.loc + t.len / 2))
                                        .attr("cy", b)
                                        .attr("r", 3));
                                const A = t
                                    .selectAll("g.orf-group")
                                    .data(e)
                                    .enter()
                                    .append("g")
                                    .attr("class", "orf-group");
                                (A.append("polygon")
                                    .attr("points", (t) => m(t))
                                    .attr("class", "svgene-orf-bg")
                                    .attr("id", (t) => `u${i}-region${f}-${x(t.locus_tag)}-svgeneorf-bg`)
                                    .style("fill", "white"),
                                    A.append("polygon")
                                        .attr("points", (t) => m(t))
                                        .attr("class", (t) => `svgene-type-${t.type} svgene-orf ${g}`)
                                        .attr("id", (t) => `u${i}-region${f}-${x(t.locus_tag)}-svgeneorf`)
                                        .attr("opacity", "1"));
                                const _ = z + d + a + 2;
                                t.selectAll("rect.svgene-resist")
                                    .data(e.filter((t) => t.resistance))
                                    .enter()
                                    .append("rect")
                                    .attr("class", "svgene-resistance")
                                    .attr("width", (t) => H(t.end) - H(t.start))
                                    .attr("height", 7)
                                    .attr("x", (t) => H(t.start))
                                    .attr("y", _);
                                for (const t of e.filter((t) => t.resistance)) {
                                    const e = `#u${i}-region${f}-${x(t.locus_tag)}-svgeneorf`;
                                    (0, s.select)(e).classed("svgene-resistance-orf", !0);
                                }
                                t.selectAll("polyline.svgene-tta-codon")
                                    .data(r.ttaCodons)
                                    .enter()
                                    .append("polyline")
                                    .attr("points", (t) => v(t, a, z, l))
                                    .attr("class", "svgene-tta-codon");
                                const C = t
                                    .selectAll("text.svgene-locustag")
                                    .data(e)
                                    .enter()
                                    .append("g")
                                    .attr("class", "svgene-locustag")
                                    .attr("id", (t) => `u${i}-region${f}-${x(t.locus_tag)}-label`);
                                (C.append("rect")
                                    .attr("x", (t) => (H(t.start) < o / 2 ? H(t.start) : H(t.end)))
                                    .attr("y", z)
                                    .attr("width", (t) => 10 * t.locus_tag.length)
                                    .attr("height", d)
                                    .attr("class", "svgene-locustag-background"),
                                    C.append("text")
                                        .attr("x", (t) => (H(t.start) < o / 2 ? H(t.start) : H(t.end)))
                                        .attr("text-anchor", (t) => (H(t.start) < o / 2 ? "start" : "end"))
                                        .attr("y", z + d)
                                        .text((t) => t.locus_tag));
                                for (const t of r.ttaCodons)
                                    for (const e of t.containedBy) {
                                        const t = `#u${i}-region${f}-${x(e)}-svgeneorf`;
                                        (0, s.select)(t).classed("contains-tta-codon", !0);
                                    }
                            })(M, k, R, N, U, n, C, F),
                            void 0 !== c.label &&
                                M.append("text")
                                    .text(c.label)
                                    .attr("class", "svgene-regionlabel")
                                    .attr("x", function () {
                                        const t = this.getComputedTextLength();
                                        return C - t - 5;
                                    })
                                    .attr("y", d)
                                    .attr("font-size", d),
                            (B = (0, i.axisBottom)(H)),
                            M.append("g")
                                .attr("class", "svgene-axis")
                                .attr("transform", `translate(0,${w - 40 - 20})`)
                                .call(B),
                            (function (t, e, n, r, i, l, c) {
                                if (!q) return;
                                O = (0, o.scaleLinear)()
                                    .domain([e, n])
                                    .range([Math.floor(0.25 * c), Math.floor(0.75 * c)])
                                    .nice();
                                (t.append("g").attr("class", "svgene-minimap"),
                                    t
                                        .append("line")
                                        .attr("x1", O(e))
                                        .attr("y1", r)
                                        .attr("x2", O(n))
                                        .attr("y2", r)
                                        .attr("class", "centerline"),
                                    t
                                        .selectAll("rect.svgene-minimap-orf")
                                        .data(l)
                                        .enter()
                                        .append("rect")
                                        .attr("width", (t) => O(t.end) - O(t.start))
                                        .attr("height", 10)
                                        .attr("x", (t) => O(t.start))
                                        .attr("y", r - 5)
                                        .attr("class", (t) => `svgene-type-${t.type} ${g} svgene-minimap-orf`)
                                        .attr("id", (t) => `svgene-minimap-orf-${x(t.locus_tag)}`));
                                const u = O(e),
                                    h = O(n),
                                    d = O(Math.min(n, e + 1e3)),
                                    f = O(Math.max(e, n - 1e3));
                                (t
                                    .append("rect")
                                    .attr("x", O(e))
                                    .attr("y", r - 20)
                                    .attr("width", O(n) - O(e))
                                    .attr("height", i)
                                    .style("cursor", "grab")
                                    .attr("id", `svgene-minimap-window-${q.anchor}`)
                                    .attr("opacity", "0.2")
                                    .style("fill", "blue")
                                    .call(
                                        (0, a.drag)().on("drag", function () {
                                            if (null === q) return;
                                            Math.max(0, s.event.x - parseInt((0, s.select)(this).attr("x"), 10));
                                            let t = parseInt((0, s.select)(this).attr("width"), 10);
                                            const e = parseInt((0, s.select)(this).attr("x"), 10) + s.event.dx;
                                            let n = y(u, e, f);
                                            (e < u
                                                ? (t = Math.max(d - u, t - (u - e)))
                                                : e + t > h && ((n = Math.min(e, f)), (t = h - n)),
                                                (n = y(u, n, f)));
                                            const r = Math.max(d, n + t);
                                            ((0, s.select)(this).attr("x", n).attr("width", t),
                                                (0, s.select)(`#svgene-minimap-grabber-a-${q.anchor}`).attr("x", n),
                                                (0, s.select)(`#svgene-minimap-grabber-b-${q.anchor}`).attr("x", r),
                                                I(O.invert(n), O.invert(r), !0));
                                        }),
                                    ),
                                    t
                                        .append("rect")
                                        .attr("x", O(e))
                                        .attr("y", r - 20)
                                        .attr("id", `svgene-minimap-grabber-a-${q.anchor}`)
                                        .attr("width", 2)
                                        .attr("height", i)
                                        .style("cursor", "ew-resize")
                                        .call(
                                            (0, a.drag)().on("drag", function (t, e) {
                                                if (!q) return;
                                                const n = parseInt(
                                                        (0, s.select)(`#svgene-minimap-grabber-b-${q.anchor}`).attr(
                                                            "x",
                                                        ),
                                                        10,
                                                    ),
                                                    r = Math.min(O(O.invert(n) - 1e3), h, Math.max(s.event.x, u));
                                                ((0, s.select)(this).transition().duration(10).attr("x", r),
                                                    (0, s.select)(`#svgene-minimap-window-${q.anchor}`)
                                                        .transition()
                                                        .duration(10)
                                                        .attr("x", r)
                                                        .attr("width", n - r),
                                                    I(O.invert(r), O.invert(n), !0));
                                            }),
                                        ),
                                    t
                                        .append("rect")
                                        .attr("x", O(n))
                                        .attr("y", r - 20)
                                        .attr("id", `svgene-minimap-grabber-b-${q.anchor}`)
                                        .attr("width", 2)
                                        .attr("height", i)
                                        .style("cursor", "ew-resize")
                                        .call(
                                            (0, a.drag)().on("drag", function (t, e) {
                                                if (!q) return;
                                                const n = parseInt(
                                                        (0, s.select)(`#svgene-minimap-grabber-a-${q.anchor}`).attr(
                                                            "x",
                                                        ),
                                                        10,
                                                    ),
                                                    r = Math.max(O(O.invert(n) + 1e3), d, Math.min(s.event.x, h));
                                                ((0, s.select)(this).transition().duration(10).attr("x", r),
                                                    (0, s.select)(`#svgene-minimap-window-${q.anchor}`)
                                                        .transition()
                                                        .duration(10)
                                                        .attr("width", r - n),
                                                    I(O.invert(n), O.invert(r), !0));
                                            }),
                                        ));
                            })(M, c.start, c.end, w - 20, 40, k, C),
                            $(".svgene-orf")
                                .mouseover(function (t) {
                                    let e = $(this).attr("id");
                                    void 0 !== e && ((e = e.replace("-svgeneorf", "-label")), $("#" + e).show());
                                })
                                .mouseout(function (t) {
                                    let e = $(this).attr("id");
                                    void 0 !== e && ((e = e.replace("-svgeneorf", "-label")), $("#" + e).hide());
                                })
                                .click(A),
                            $(".svgene-textarea").click((t) => t.stopPropagation()),
                            $(".legend-selector").unbind("click").click(P),
                            $(".zoom-in").unbind("click").click(D),
                            $(".zoom-reset")
                                .unbind("click")
                                .click(() => E()),
                            1 === j &&
                                $(document).keyup((t) => {
                                    const e = t.keyCode;
                                    82 === e ? E() : 90 === e && D();
                                }));
                    }),
                    (e.tag_to_id = x),
                    (e.locusSelectionHandler = function (t) {
                        const e = this.getAttribute("data-locus");
                        e && f([e], t.ctrlKey || t.metaKey);
                    }),
                    (e.zoom_to_selection = D),
                    (e.resetView = E),
                    (e.resetZoom = function () {
                        null !== q && I(q.start, q.end);
                    }));
                let B = null,
                    O = (0, o.scaleLinear)().domain([0, 100]).range([0, 100]),
                    F = 0,
                    z = 0,
                    q = null,
                    H = (0, o.scaleLinear)().domain([0, 100]).range([0, 100]),
                    j = 0;
            },
            130: (t, e) => {
                "use strict";
                (Object.defineProperty(e, "__esModule", { value: !0 }),
                    (e.replaceWildcards = e.replaceWildcardsInText = void 0));
                const n = "wildcard-container",
                    r = /@![^!]*!@/g;
                function i(t, e) {
                    for (const n of t.match(r) || []) {
                        const r = `${e[n.substring(2, n.length - 2)] || ""}`;
                        t = t.replace(n, r);
                    }
                    return t;
                }
                ((e.replaceWildcardsInText = i),
                    (e.replaceWildcards = function (t, e) {
                        $(t)
                            .find(`.${n}`)
                            .each(function () {
                                this.classList.remove(n);
                                const t = this.getAttribute("data-locus");
                                if (void 0 === t) return;
                                const r = e.orfs.filter((e) => e.locus_tag === t)[0];
                                for (const t of (this.getAttribute("data-wildcard-attrs") || "").split(" ")) {
                                    const e = this.getAttribute(t);
                                    e && this.setAttribute(t, i(e, r));
                                }
                            });
                    }));
            },
            216: (t, e, n) => {
                "use strict";
                (n.r(e), n.d(e, { axisBottom: () => y, axisLeft: () => b, axisRight: () => v, axisTop: () => m }));
                var r = Array.prototype.slice;
                function i(t) {
                    return t;
                }
                var a = 1,
                    o = 2,
                    s = 3,
                    l = 4,
                    c = 1e-6;
                function u(t) {
                    return "translate(" + (t + 0.5) + ",0)";
                }
                function h(t) {
                    return "translate(0," + (t + 0.5) + ")";
                }
                function d(t) {
                    return function (e) {
                        return +t(e);
                    };
                }
                function g(t) {
                    var e = Math.max(0, t.bandwidth() - 1) / 2;
                    return (
                        t.round() && (e = Math.round(e)),
                        function (n) {
                            return +t(n) + e;
                        }
                    );
                }
                function f() {
                    return !this.__axis;
                }
                function p(t, e) {
                    var n = [],
                        p = null,
                        m = null,
                        v = 6,
                        y = 6,
                        b = 3,
                        x = t === a || t === l ? -1 : 1,
                        w = t === l || t === o ? "x" : "y",
                        $ = t === a || t === s ? u : h;
                    function A(r) {
                        var u = null == p ? (e.ticks ? e.ticks.apply(e, n) : e.domain()) : p,
                            h = null == m ? (e.tickFormat ? e.tickFormat.apply(e, n) : i) : m,
                            A = Math.max(v, 0) + b,
                            _ = e.range(),
                            C = +_[0] + 0.5,
                            S = +_[_.length - 1] + 0.5,
                            M = (e.bandwidth ? g : d)(e.copy()),
                            T = r.selection ? r.selection() : r,
                            k = T.selectAll(".domain").data([null]),
                            R = T.selectAll(".tick").data(u, e).order(),
                            N = R.exit(),
                            L = R.enter().append("g").attr("class", "tick"),
                            P = R.select("line"),
                            D = R.select("text");
                        ((k = k.merge(
                            k.enter().insert("path", ".tick").attr("class", "domain").attr("stroke", "currentColor"),
                        )),
                            (R = R.merge(L)),
                            (P = P.merge(
                                L.append("line")
                                    .attr("stroke", "currentColor")
                                    .attr(w + "2", x * v),
                            )),
                            (D = D.merge(
                                L.append("text")
                                    .attr("fill", "currentColor")
                                    .attr(w, x * A)
                                    .attr("dy", t === a ? "0em" : t === s ? "0.71em" : "0.32em"),
                            )),
                            r !== T &&
                                ((k = k.transition(r)),
                                (R = R.transition(r)),
                                (P = P.transition(r)),
                                (D = D.transition(r)),
                                (N = N.transition(r)
                                    .attr("opacity", c)
                                    .attr("transform", function (t) {
                                        return isFinite((t = M(t))) ? $(t) : this.getAttribute("transform");
                                    })),
                                L.attr("opacity", c).attr("transform", function (t) {
                                    var e = this.parentNode.__axis;
                                    return $(e && isFinite((e = e(t))) ? e : M(t));
                                })),
                            N.remove(),
                            k.attr(
                                "d",
                                t === l || t == o
                                    ? y
                                        ? "M" + x * y + "," + C + "H0.5V" + S + "H" + x * y
                                        : "M0.5," + C + "V" + S
                                    : y
                                      ? "M" + C + "," + x * y + "V0.5H" + S + "V" + x * y
                                      : "M" + C + ",0.5H" + S,
                            ),
                            R.attr("opacity", 1).attr("transform", function (t) {
                                return $(M(t));
                            }),
                            P.attr(w + "2", x * v),
                            D.attr(w, x * A).text(h),
                            T.filter(f)
                                .attr("fill", "none")
                                .attr("font-size", 10)
                                .attr("font-family", "sans-serif")
                                .attr("text-anchor", t === o ? "start" : t === l ? "end" : "middle"),
                            T.each(function () {
                                this.__axis = M;
                            }));
                    }
                    return (
                        (A.scale = function (t) {
                            return arguments.length ? ((e = t), A) : e;
                        }),
                        (A.ticks = function () {
                            return ((n = r.call(arguments)), A);
                        }),
                        (A.tickArguments = function (t) {
                            return arguments.length ? ((n = null == t ? [] : r.call(t)), A) : n.slice();
                        }),
                        (A.tickValues = function (t) {
                            return arguments.length ? ((p = null == t ? null : r.call(t)), A) : p && p.slice();
                        }),
                        (A.tickFormat = function (t) {
                            return arguments.length ? ((m = t), A) : m;
                        }),
                        (A.tickSize = function (t) {
                            return arguments.length ? ((v = y = +t), A) : v;
                        }),
                        (A.tickSizeInner = function (t) {
                            return arguments.length ? ((v = +t), A) : v;
                        }),
                        (A.tickSizeOuter = function (t) {
                            return arguments.length ? ((y = +t), A) : y;
                        }),
                        (A.tickPadding = function (t) {
                            return arguments.length ? ((b = +t), A) : b;
                        }),
                        A
                    );
                }
                function m(t) {
                    return p(a, t);
                }
                function v(t) {
                    return p(o, t);
                }
                function y(t) {
                    return p(s, t);
                }
                function b(t) {
                    return p(l, t);
                }
            },
            650: (t, e, n) => {
                "use strict";
                function r(t, e, n) {
                    ((t.prototype = e.prototype = n), (n.constructor = t));
                }
                function i(t, e) {
                    var n = Object.create(t.prototype);
                    for (var r in e) n[r] = e[r];
                    return n;
                }
                function a() {}
                n.d(e, { ZP: () => w, B8: () => _ });
                var o = 0.7,
                    s = 1 / o,
                    l = "\\s*([+-]?\\d+)\\s*",
                    c = "\\s*([+-]?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)\\s*",
                    u = "\\s*([+-]?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?)%\\s*",
                    h = /^#([0-9a-f]{3,8})$/,
                    d = new RegExp("^rgb\\(" + [l, l, l] + "\\)$"),
                    g = new RegExp("^rgb\\(" + [u, u, u] + "\\)$"),
                    f = new RegExp("^rgba\\(" + [l, l, l, c] + "\\)$"),
                    p = new RegExp("^rgba\\(" + [u, u, u, c] + "\\)$"),
                    m = new RegExp("^hsl\\(" + [c, u, u] + "\\)$"),
                    v = new RegExp("^hsla\\(" + [c, u, u, c] + "\\)$"),
                    y = {
                        aliceblue: 15792383,
                        antiquewhite: 16444375,
                        aqua: 65535,
                        aquamarine: 8388564,
                        azure: 15794175,
                        beige: 16119260,
                        bisque: 16770244,
                        black: 0,
                        blanchedalmond: 16772045,
                        blue: 255,
                        blueviolet: 9055202,
                        brown: 10824234,
                        burlywood: 14596231,
                        cadetblue: 6266528,
                        chartreuse: 8388352,
                        chocolate: 13789470,
                        coral: 16744272,
                        cornflowerblue: 6591981,
                        cornsilk: 16775388,
                        crimson: 14423100,
                        cyan: 65535,
                        darkblue: 139,
                        darkcyan: 35723,
                        darkgoldenrod: 12092939,
                        darkgray: 11119017,
                        darkgreen: 25600,
                        darkgrey: 11119017,
                        darkkhaki: 12433259,
                        darkmagenta: 9109643,
                        darkolivegreen: 5597999,
                        darkorange: 16747520,
                        darkorchid: 10040012,
                        darkred: 9109504,
                        darksalmon: 15308410,
                        darkseagreen: 9419919,
                        darkslateblue: 4734347,
                        darkslategray: 3100495,
                        darkslategrey: 3100495,
                        darkturquoise: 52945,
                        darkviolet: 9699539,
                        deeppink: 16716947,
                        deepskyblue: 49151,
                        dimgray: 6908265,
                        dimgrey: 6908265,
                        dodgerblue: 2003199,
                        firebrick: 11674146,
                        floralwhite: 16775920,
                        forestgreen: 2263842,
                        fuchsia: 16711935,
                        gainsboro: 14474460,
                        ghostwhite: 16316671,
                        gold: 16766720,
                        goldenrod: 14329120,
                        gray: 8421504,
                        green: 32768,
                        greenyellow: 11403055,
                        grey: 8421504,
                        honeydew: 15794160,
                        hotpink: 16738740,
                        indianred: 13458524,
                        indigo: 4915330,
                        ivory: 16777200,
                        khaki: 15787660,
                        lavender: 15132410,
                        lavenderblush: 16773365,
                        lawngreen: 8190976,
                        lemonchiffon: 16775885,
                        lightblue: 11393254,
                        lightcoral: 15761536,
                        lightcyan: 14745599,
                        lightgoldenrodyellow: 16448210,
                        lightgray: 13882323,
                        lightgreen: 9498256,
                        lightgrey: 13882323,
                        lightpink: 16758465,
                        lightsalmon: 16752762,
                        lightseagreen: 2142890,
                        lightskyblue: 8900346,
                        lightslategray: 7833753,
                        lightslategrey: 7833753,
                        lightsteelblue: 11584734,
                        lightyellow: 16777184,
                        lime: 65280,
                        limegreen: 3329330,
                        linen: 16445670,
                        magenta: 16711935,
                        maroon: 8388608,
                        mediumaquamarine: 6737322,
                        mediumblue: 205,
                        mediumorchid: 12211667,
                        mediumpurple: 9662683,
                        mediumseagreen: 3978097,
                        mediumslateblue: 8087790,
                        mediumspringgreen: 64154,
                        mediumturquoise: 4772300,
                        mediumvioletred: 13047173,
                        midnightblue: 1644912,
                        mintcream: 16121850,
                        mistyrose: 16770273,
                        moccasin: 16770229,
                        navajowhite: 16768685,
                        navy: 128,
                        oldlace: 16643558,
                        olive: 8421376,
                        olivedrab: 7048739,
                        orange: 16753920,
                        orangered: 16729344,
                        orchid: 14315734,
                        palegoldenrod: 15657130,
                        palegreen: 10025880,
                        paleturquoise: 11529966,
                        palevioletred: 14381203,
                        papayawhip: 16773077,
                        peachpuff: 16767673,
                        peru: 13468991,
                        pink: 16761035,
                        plum: 14524637,
                        powderblue: 11591910,
                        purple: 8388736,
                        rebeccapurple: 6697881,
                        red: 16711680,
                        rosybrown: 12357519,
                        royalblue: 4286945,
                        saddlebrown: 9127187,
                        salmon: 16416882,
                        sandybrown: 16032864,
                        seagreen: 3050327,
                        seashell: 16774638,
                        sienna: 10506797,
                        silver: 12632256,
                        skyblue: 8900331,
                        slateblue: 6970061,
                        slategray: 7372944,
                        slategrey: 7372944,
                        snow: 16775930,
                        springgreen: 65407,
                        steelblue: 4620980,
                        tan: 13808780,
                        teal: 32896,
                        thistle: 14204888,
                        tomato: 16737095,
                        turquoise: 4251856,
                        violet: 15631086,
                        wheat: 16113331,
                        white: 16777215,
                        whitesmoke: 16119285,
                        yellow: 16776960,
                        yellowgreen: 10145074,
                    };
                function b() {
                    return this.rgb().formatHex();
                }
                function x() {
                    return this.rgb().formatRgb();
                }
                function w(t) {
                    var e, n;
                    return (
                        (t = (t + "").trim().toLowerCase()),
                        (e = h.exec(t))
                            ? ((n = e[1].length),
                              (e = parseInt(e[1], 16)),
                              6 === n
                                  ? $(e)
                                  : 3 === n
                                    ? new C(
                                          ((e >> 8) & 15) | ((e >> 4) & 240),
                                          ((e >> 4) & 15) | (240 & e),
                                          ((15 & e) << 4) | (15 & e),
                                          1,
                                      )
                                    : 8 === n
                                      ? A((e >> 24) & 255, (e >> 16) & 255, (e >> 8) & 255, (255 & e) / 255)
                                      : 4 === n
                                        ? A(
                                              ((e >> 12) & 15) | ((e >> 8) & 240),
                                              ((e >> 8) & 15) | ((e >> 4) & 240),
                                              ((e >> 4) & 15) | (240 & e),
                                              (((15 & e) << 4) | (15 & e)) / 255,
                                          )
                                        : null)
                            : (e = d.exec(t))
                              ? new C(e[1], e[2], e[3], 1)
                              : (e = g.exec(t))
                                ? new C((255 * e[1]) / 100, (255 * e[2]) / 100, (255 * e[3]) / 100, 1)
                                : (e = f.exec(t))
                                  ? A(e[1], e[2], e[3], e[4])
                                  : (e = p.exec(t))
                                    ? A((255 * e[1]) / 100, (255 * e[2]) / 100, (255 * e[3]) / 100, e[4])
                                    : (e = m.exec(t))
                                      ? k(e[1], e[2] / 100, e[3] / 100, 1)
                                      : (e = v.exec(t))
                                        ? k(e[1], e[2] / 100, e[3] / 100, e[4])
                                        : y.hasOwnProperty(t)
                                          ? $(y[t])
                                          : "transparent" === t
                                            ? new C(NaN, NaN, NaN, 0)
                                            : null
                    );
                }
                function $(t) {
                    return new C((t >> 16) & 255, (t >> 8) & 255, 255 & t, 1);
                }
                function A(t, e, n, r) {
                    return (r <= 0 && (t = e = n = NaN), new C(t, e, n, r));
                }
                function _(t, e, n, r) {
                    return 1 === arguments.length
                        ? ((i = t) instanceof a || (i = w(i)),
                          i ? new C((i = i.rgb()).r, i.g, i.b, i.opacity) : new C())
                        : new C(t, e, n, null == r ? 1 : r);
                    var i;
                }
                function C(t, e, n, r) {
                    ((this.r = +t), (this.g = +e), (this.b = +n), (this.opacity = +r));
                }
                function S() {
                    return "#" + T(this.r) + T(this.g) + T(this.b);
                }
                function M() {
                    var t = this.opacity;
                    return (
                        (1 === (t = isNaN(t) ? 1 : Math.max(0, Math.min(1, t))) ? "rgb(" : "rgba(") +
                        Math.max(0, Math.min(255, Math.round(this.r) || 0)) +
                        ", " +
                        Math.max(0, Math.min(255, Math.round(this.g) || 0)) +
                        ", " +
                        Math.max(0, Math.min(255, Math.round(this.b) || 0)) +
                        (1 === t ? ")" : ", " + t + ")")
                    );
                }
                function T(t) {
                    return ((t = Math.max(0, Math.min(255, Math.round(t) || 0))) < 16 ? "0" : "") + t.toString(16);
                }
                function k(t, e, n, r) {
                    return (
                        r <= 0 ? (t = e = n = NaN) : n <= 0 || n >= 1 ? (t = e = NaN) : e <= 0 && (t = NaN),
                        new N(t, e, n, r)
                    );
                }
                function R(t) {
                    if (t instanceof N) return new N(t.h, t.s, t.l, t.opacity);
                    if ((t instanceof a || (t = w(t)), !t)) return new N();
                    if (t instanceof N) return t;
                    var e = (t = t.rgb()).r / 255,
                        n = t.g / 255,
                        r = t.b / 255,
                        i = Math.min(e, n, r),
                        o = Math.max(e, n, r),
                        s = NaN,
                        l = o - i,
                        c = (o + i) / 2;
                    return (
                        l
                            ? ((s = e === o ? (n - r) / l + 6 * (n < r) : n === o ? (r - e) / l + 2 : (e - n) / l + 4),
                              (l /= c < 0.5 ? o + i : 2 - o - i),
                              (s *= 60))
                            : (l = c > 0 && c < 1 ? 0 : s),
                        new N(s, l, c, t.opacity)
                    );
                }
                function N(t, e, n, r) {
                    ((this.h = +t), (this.s = +e), (this.l = +n), (this.opacity = +r));
                }
                function L(t, e, n) {
                    return (
                        255 *
                        (t < 60 ? e + ((n - e) * t) / 60 : t < 180 ? n : t < 240 ? e + ((n - e) * (240 - t)) / 60 : e)
                    );
                }
                (r(a, w, {
                    copy: function (t) {
                        return Object.assign(new this.constructor(), this, t);
                    },
                    displayable: function () {
                        return this.rgb().displayable();
                    },
                    hex: b,
                    formatHex: b,
                    formatHsl: function () {
                        return R(this).formatHsl();
                    },
                    formatRgb: x,
                    toString: x,
                }),
                    r(
                        C,
                        _,
                        i(a, {
                            brighter: function (t) {
                                return (
                                    (t = null == t ? s : Math.pow(s, t)),
                                    new C(this.r * t, this.g * t, this.b * t, this.opacity)
                                );
                            },
                            darker: function (t) {
                                return (
                                    (t = null == t ? o : Math.pow(o, t)),
                                    new C(this.r * t, this.g * t, this.b * t, this.opacity)
                                );
                            },
                            rgb: function () {
                                return this;
                            },
                            displayable: function () {
                                return (
                                    -0.5 <= this.r &&
                                    this.r < 255.5 &&
                                    -0.5 <= this.g &&
                                    this.g < 255.5 &&
                                    -0.5 <= this.b &&
                                    this.b < 255.5 &&
                                    0 <= this.opacity &&
                                    this.opacity <= 1
                                );
                            },
                            hex: S,
                            formatHex: S,
                            formatRgb: M,
                            toString: M,
                        }),
                    ),
                    r(
                        N,
                        function (t, e, n, r) {
                            return 1 === arguments.length ? R(t) : new N(t, e, n, null == r ? 1 : r);
                        },
                        i(a, {
                            brighter: function (t) {
                                return (
                                    (t = null == t ? s : Math.pow(s, t)),
                                    new N(this.h, this.s, this.l * t, this.opacity)
                                );
                            },
                            darker: function (t) {
                                return (
                                    (t = null == t ? o : Math.pow(o, t)),
                                    new N(this.h, this.s, this.l * t, this.opacity)
                                );
                            },
                            rgb: function () {
                                var t = (this.h % 360) + 360 * (this.h < 0),
                                    e = isNaN(t) || isNaN(this.s) ? 0 : this.s,
                                    n = this.l,
                                    r = n + (n < 0.5 ? n : 1 - n) * e,
                                    i = 2 * n - r;
                                return new C(
                                    L(t >= 240 ? t - 240 : t + 120, i, r),
                                    L(t, i, r),
                                    L(t < 120 ? t + 240 : t - 120, i, r),
                                    this.opacity,
                                );
                            },
                            displayable: function () {
                                return (
                                    ((0 <= this.s && this.s <= 1) || isNaN(this.s)) &&
                                    0 <= this.l &&
                                    this.l <= 1 &&
                                    0 <= this.opacity &&
                                    this.opacity <= 1
                                );
                            },
                            formatHsl: function () {
                                var t = this.opacity;
                                return (
                                    (1 === (t = isNaN(t) ? 1 : Math.max(0, Math.min(1, t))) ? "hsl(" : "hsla(") +
                                    (this.h || 0) +
                                    ", " +
                                    100 * (this.s || 0) +
                                    "%, " +
                                    100 * (this.l || 0) +
                                    "%" +
                                    (1 === t ? ")" : ", " + t + ")")
                                );
                            },
                        }),
                    ));
            },
            594: (t, e, n) => {
                "use strict";
                n.d(e, { W: () => l });
                var r = { value: function () {} };
                function i() {
                    for (var t, e = 0, n = arguments.length, r = {}; e < n; ++e) {
                        if (!(t = arguments[e] + "") || t in r) throw new Error("illegal type: " + t);
                        r[t] = [];
                    }
                    return new a(r);
                }
                function a(t) {
                    this._ = t;
                }
                function o(t, e) {
                    for (var n, r = 0, i = t.length; r < i; ++r) if ((n = t[r]).name === e) return n.value;
                }
                function s(t, e, n) {
                    for (var i = 0, a = t.length; i < a; ++i)
                        if (t[i].name === e) {
                            ((t[i] = r), (t = t.slice(0, i).concat(t.slice(i + 1))));
                            break;
                        }
                    return (null != n && t.push({ name: e, value: n }), t);
                }
                a.prototype = i.prototype = {
                    constructor: a,
                    on: function (t, e) {
                        var n,
                            r,
                            i = this._,
                            a =
                                ((r = i),
                                (t + "")
                                    .trim()
                                    .split(/^|\s+/)
                                    .map(function (t) {
                                        var e = "",
                                            n = t.indexOf(".");
                                        if (
                                            (n >= 0 && ((e = t.slice(n + 1)), (t = t.slice(0, n))),
                                            t && !r.hasOwnProperty(t))
                                        )
                                            throw new Error("unknown type: " + t);
                                        return { type: t, name: e };
                                    })),
                            l = -1,
                            c = a.length;
                        if (!(arguments.length < 2)) {
                            if (null != e && "function" != typeof e) throw new Error("invalid callback: " + e);
                            for (; ++l < c; )
                                if ((n = (t = a[l]).type)) i[n] = s(i[n], t.name, e);
                                else if (null == e) for (n in i) i[n] = s(i[n], t.name, null);
                            return this;
                        }
                        for (; ++l < c; ) if ((n = (t = a[l]).type) && (n = o(i[n], t.name))) return n;
                    },
                    copy: function () {
                        var t = {},
                            e = this._;
                        for (var n in e) t[n] = e[n].slice();
                        return new a(t);
                    },
                    call: function (t, e) {
                        if ((n = arguments.length - 2) > 0)
                            for (var n, r, i = new Array(n), a = 0; a < n; ++a) i[a] = arguments[a + 2];
                        if (!this._.hasOwnProperty(t)) throw new Error("unknown type: " + t);
                        for (a = 0, n = (r = this._[t]).length; a < n; ++a) r[a].value.apply(e, i);
                    },
                    apply: function (t, e, n) {
                        if (!this._.hasOwnProperty(t)) throw new Error("unknown type: " + t);
                        for (var r = this._[t], i = 0, a = r.length; i < a; ++i) r[i].value.apply(e, n);
                    },
                };
                const l = i;
            },
            305: (t, e, n) => {
                "use strict";
                (n.r(e), n.d(e, { drag: () => p, dragDisable: () => s, dragEnable: () => l }));
                var r = n(594),
                    i = n(677);
                function a() {
                    i.event.stopImmediatePropagation();
                }
                function o() {
                    (i.event.preventDefault(), i.event.stopImmediatePropagation());
                }
                function s(t) {
                    var e = t.document.documentElement,
                        n = (0, i.select)(t).on("dragstart.drag", o, !0);
                    "onselectstart" in e
                        ? n.on("selectstart.drag", o, !0)
                        : ((e.__noselect = e.style.MozUserSelect), (e.style.MozUserSelect = "none"));
                }
                function l(t, e) {
                    var n = t.document.documentElement,
                        r = (0, i.select)(t).on("dragstart.drag", null);
                    (e &&
                        (r.on("click.drag", o, !0),
                        setTimeout(function () {
                            r.on("click.drag", null);
                        }, 0)),
                        "onselectstart" in n
                            ? r.on("selectstart.drag", null)
                            : ((n.style.MozUserSelect = n.__noselect), delete n.__noselect));
                }
                function c(t) {
                    return function () {
                        return t;
                    };
                }
                function u(t, e, n, r, i, a, o, s, l, c) {
                    ((this.target = t),
                        (this.type = e),
                        (this.subject = n),
                        (this.identifier = r),
                        (this.active = i),
                        (this.x = a),
                        (this.y = o),
                        (this.dx = s),
                        (this.dy = l),
                        (this._ = c));
                }
                function h() {
                    return !i.event.button;
                }
                function d() {
                    return this.parentNode;
                }
                function g(t) {
                    return null == t ? { x: i.event.x, y: i.event.y } : t;
                }
                function f() {
                    return "ontouchstart" in this;
                }
                function p() {
                    var t,
                        e,
                        n,
                        p,
                        m = h,
                        v = d,
                        y = g,
                        b = f,
                        x = {},
                        w = (0, r.W)("start", "drag", "end"),
                        $ = 0,
                        A = 0;
                    function _(t) {
                        t.on("mousedown.drag", C)
                            .filter(b)
                            .on("touchstart.drag", T)
                            .on("touchmove.drag", k)
                            .on("touchend.drag touchcancel.drag", R)
                            .style("touch-action", "none")
                            .style("-webkit-tap-highlight-color", "rgba(0,0,0,0)");
                    }
                    function C() {
                        if (!p && m.apply(this, arguments)) {
                            var r = N("mouse", v.apply(this, arguments), i.mouse, this, arguments);
                            r &&
                                ((0, i.select)(i.event.view).on("mousemove.drag", S, !0).on("mouseup.drag", M, !0),
                                s(i.event.view),
                                a(),
                                (n = !1),
                                (t = i.event.clientX),
                                (e = i.event.clientY),
                                r("start"));
                        }
                    }
                    function S() {
                        if ((o(), !n)) {
                            var r = i.event.clientX - t,
                                a = i.event.clientY - e;
                            n = r * r + a * a > A;
                        }
                        x.mouse("drag");
                    }
                    function M() {
                        ((0, i.select)(i.event.view).on("mousemove.drag mouseup.drag", null),
                            l(i.event.view, n),
                            o(),
                            x.mouse("end"));
                    }
                    function T() {
                        if (m.apply(this, arguments)) {
                            var t,
                                e,
                                n = i.event.changedTouches,
                                r = v.apply(this, arguments),
                                o = n.length;
                            for (t = 0; t < o; ++t)
                                (e = N(n[t].identifier, r, i.touch, this, arguments)) && (a(), e("start"));
                        }
                    }
                    function k() {
                        var t,
                            e,
                            n = i.event.changedTouches,
                            r = n.length;
                        for (t = 0; t < r; ++t) (e = x[n[t].identifier]) && (o(), e("drag"));
                    }
                    function R() {
                        var t,
                            e,
                            n = i.event.changedTouches,
                            r = n.length;
                        for (
                            p && clearTimeout(p),
                                p = setTimeout(function () {
                                    p = null;
                                }, 500),
                                t = 0;
                            t < r;
                            ++t
                        )
                            (e = x[n[t].identifier]) && (a(), e("end"));
                    }
                    function N(t, e, n, r, a) {
                        var o,
                            s,
                            l,
                            c = n(e, t),
                            h = w.copy();
                        if (
                            (0, i.customEvent)(new u(_, "beforestart", o, t, $, c[0], c[1], 0, 0, h), function () {
                                return (
                                    null != (i.event.subject = o = y.apply(r, a)) &&
                                    ((s = o.x - c[0] || 0), (l = o.y - c[1] || 0), !0)
                                );
                            })
                        )
                            return function d(g) {
                                var f,
                                    p = c;
                                switch (g) {
                                    case "start":
                                        ((x[t] = d), (f = $++));
                                        break;
                                    case "end":
                                        (delete x[t], --$);
                                    case "drag":
                                        ((c = n(e, t)), (f = $));
                                }
                                (0, i.customEvent)(
                                    new u(_, g, o, t, f, c[0] + s, c[1] + l, c[0] - p[0], c[1] - p[1], h),
                                    h.apply,
                                    h,
                                    [g, r, a],
                                );
                            };
                    }
                    return (
                        (_.filter = function (t) {
                            return arguments.length ? ((m = "function" == typeof t ? t : c(!!t)), _) : m;
                        }),
                        (_.container = function (t) {
                            return arguments.length ? ((v = "function" == typeof t ? t : c(t)), _) : v;
                        }),
                        (_.subject = function (t) {
                            return arguments.length ? ((y = "function" == typeof t ? t : c(t)), _) : y;
                        }),
                        (_.touchable = function (t) {
                            return arguments.length ? ((b = "function" == typeof t ? t : c(!!t)), _) : b;
                        }),
                        (_.on = function () {
                            var t = w.on.apply(w, arguments);
                            return t === w ? _ : t;
                        }),
                        (_.clickDistance = function (t) {
                            return arguments.length ? ((A = (t = +t) * t), _) : Math.sqrt(A);
                        }),
                        _
                    );
                }
                u.prototype.on = function () {
                    var t = this._.on.apply(this._, arguments);
                    return t === this._ ? this : t;
                };
            },
            302: (t, e, n) => {
                "use strict";
                function r(t) {
                    return function () {
                        return t;
                    };
                }
                n.d(e, { Z: () => r });
            },
            626: (t, e, n) => {
                "use strict";
                function r(t, e) {
                    return (
                        (t = +t),
                        (e = +e),
                        function (n) {
                            return t * (1 - n) + e * n;
                        }
                    );
                }
                n.d(e, { Z: () => r });
            },
            765: (t, e, n) => {
                "use strict";
                n.d(e, { ZP: () => s });
                var r = n(650);
                function i(t, e, n, r, i) {
                    var a = t * t,
                        o = a * t;
                    return (
                        ((1 - 3 * t + 3 * a - o) * e +
                            (4 - 6 * a + 3 * o) * n +
                            (1 + 3 * t + 3 * a - 3 * o) * r +
                            o * i) /
                        6
                    );
                }
                var a = n(302);
                function o(t, e) {
                    var n = e - t;
                    return n
                        ? (function (t, e) {
                              return function (n) {
                                  return t + n * e;
                              };
                          })(t, n)
                        : (0, a.Z)(isNaN(t) ? e : t);
                }
                const s = (function t(e) {
                    var n = (function (t) {
                        return 1 == (t = +t)
                            ? o
                            : function (e, n) {
                                  return n - e
                                      ? (function (t, e, n) {
                                            return (
                                                (t = Math.pow(t, n)),
                                                (e = Math.pow(e, n) - t),
                                                (n = 1 / n),
                                                function (r) {
                                                    return Math.pow(t + r * e, n);
                                                }
                                            );
                                        })(e, n, t)
                                      : (0, a.Z)(isNaN(e) ? n : e);
                              };
                    })(e);
                    function i(t, e) {
                        var i = n((t = (0, r.B8)(t)).r, (e = (0, r.B8)(e)).r),
                            a = n(t.g, e.g),
                            s = n(t.b, e.b),
                            l = o(t.opacity, e.opacity);
                        return function (e) {
                            return ((t.r = i(e)), (t.g = a(e)), (t.b = s(e)), (t.opacity = l(e)), t + "");
                        };
                    }
                    return ((i.gamma = t), i);
                })(1);
                function l(t) {
                    return function (e) {
                        var n,
                            i,
                            a = e.length,
                            o = new Array(a),
                            s = new Array(a),
                            l = new Array(a);
                        for (n = 0; n < a; ++n)
                            ((i = (0, r.B8)(e[n])), (o[n] = i.r || 0), (s[n] = i.g || 0), (l[n] = i.b || 0));
                        return (
                            (o = t(o)),
                            (s = t(s)),
                            (l = t(l)),
                            (i.opacity = 1),
                            function (t) {
                                return ((i.r = o(t)), (i.g = s(t)), (i.b = l(t)), i + "");
                            }
                        );
                    };
                }
                (l(function (t) {
                    var e = t.length - 1;
                    return function (n) {
                        var r = n <= 0 ? (n = 0) : n >= 1 ? ((n = 1), e - 1) : Math.floor(n * e),
                            a = t[r],
                            o = t[r + 1],
                            s = r > 0 ? t[r - 1] : 2 * a - o,
                            l = r < e - 1 ? t[r + 2] : 2 * o - a;
                        return i((n - r / e) * e, s, a, o, l);
                    };
                }),
                    l(function (t) {
                        var e = t.length;
                        return function (n) {
                            var r = Math.floor(((n %= 1) < 0 ? ++n : n) * e),
                                a = t[(r + e - 1) % e],
                                o = t[r % e],
                                s = t[(r + 1) % e],
                                l = t[(r + 2) % e];
                            return i((n - r / e) * e, a, o, s, l);
                        };
                    }));
            },
            843: (t, e, n) => {
                "use strict";
                n.d(e, { Z: () => o });
                var r = n(626),
                    i = /[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g,
                    a = new RegExp(i.source, "g");
                function o(t, e) {
                    var n,
                        o,
                        s,
                        l = (i.lastIndex = a.lastIndex = 0),
                        c = -1,
                        u = [],
                        h = [];
                    for (t += "", e += ""; (n = i.exec(t)) && (o = a.exec(e)); )
                        ((s = o.index) > l && ((s = e.slice(l, s)), u[c] ? (u[c] += s) : (u[++c] = s)),
                            (n = n[0]) === (o = o[0])
                                ? u[c]
                                    ? (u[c] += o)
                                    : (u[++c] = o)
                                : ((u[++c] = null), h.push({ i: c, x: (0, r.Z)(n, o) })),
                            (l = a.lastIndex));
                    return (
                        l < e.length && ((s = e.slice(l)), u[c] ? (u[c] += s) : (u[++c] = s)),
                        u.length < 2
                            ? h[0]
                                ? (function (t) {
                                      return function (e) {
                                          return t(e) + "";
                                      };
                                  })(h[0].x)
                                : (function (t) {
                                      return function () {
                                          return t;
                                      };
                                  })(e)
                            : ((e = h.length),
                              function (t) {
                                  for (var n, r = 0; r < e; ++r) u[(n = h[r]).i] = n.x(t);
                                  return u.join("");
                              })
                    );
                }
            },
            777: (t, e, n) => {
                "use strict";
                (n.r(e), n.d(e, { path: () => c }));
                var r = Math.PI,
                    i = 2 * r,
                    a = 1e-6,
                    o = i - a;
                function s() {
                    ((this._x0 = this._y0 = this._x1 = this._y1 = null), (this._ = ""));
                }
                function l() {
                    return new s();
                }
                s.prototype = l.prototype = {
                    constructor: s,
                    moveTo: function (t, e) {
                        this._ += "M" + (this._x0 = this._x1 = +t) + "," + (this._y0 = this._y1 = +e);
                    },
                    closePath: function () {
                        null !== this._x1 && ((this._x1 = this._x0), (this._y1 = this._y0), (this._ += "Z"));
                    },
                    lineTo: function (t, e) {
                        this._ += "L" + (this._x1 = +t) + "," + (this._y1 = +e);
                    },
                    quadraticCurveTo: function (t, e, n, r) {
                        this._ += "Q" + +t + "," + +e + "," + (this._x1 = +n) + "," + (this._y1 = +r);
                    },
                    bezierCurveTo: function (t, e, n, r, i, a) {
                        this._ +=
                            "C" + +t + "," + +e + "," + +n + "," + +r + "," + (this._x1 = +i) + "," + (this._y1 = +a);
                    },
                    arcTo: function (t, e, n, i, o) {
                        ((t = +t), (e = +e), (n = +n), (i = +i), (o = +o));
                        var s = this._x1,
                            l = this._y1,
                            c = n - t,
                            u = i - e,
                            h = s - t,
                            d = l - e,
                            g = h * h + d * d;
                        if (o < 0) throw new Error("negative radius: " + o);
                        if (null === this._x1) this._ += "M" + (this._x1 = t) + "," + (this._y1 = e);
                        else if (g > a)
                            if (Math.abs(d * c - u * h) > a && o) {
                                var f = n - s,
                                    p = i - l,
                                    m = c * c + u * u,
                                    v = f * f + p * p,
                                    y = Math.sqrt(m),
                                    b = Math.sqrt(g),
                                    x = o * Math.tan((r - Math.acos((m + g - v) / (2 * y * b))) / 2),
                                    w = x / b,
                                    $ = x / y;
                                (Math.abs(w - 1) > a && (this._ += "L" + (t + w * h) + "," + (e + w * d)),
                                    (this._ +=
                                        "A" +
                                        o +
                                        "," +
                                        o +
                                        ",0,0," +
                                        +(d * f > h * p) +
                                        "," +
                                        (this._x1 = t + $ * c) +
                                        "," +
                                        (this._y1 = e + $ * u)));
                            } else this._ += "L" + (this._x1 = t) + "," + (this._y1 = e);
                    },
                    arc: function (t, e, n, s, l, c) {
                        ((t = +t), (e = +e));
                        var u = (n = +n) * Math.cos(s),
                            h = n * Math.sin(s),
                            d = t + u,
                            g = e + h,
                            f = 1 ^ c,
                            p = c ? s - l : l - s;
                        if (n < 0) throw new Error("negative radius: " + n);
                        (null === this._x1
                            ? (this._ += "M" + d + "," + g)
                            : (Math.abs(this._x1 - d) > a || Math.abs(this._y1 - g) > a) &&
                              (this._ += "L" + d + "," + g),
                            n &&
                                (p < 0 && (p = (p % i) + i),
                                p > o
                                    ? (this._ +=
                                          "A" +
                                          n +
                                          "," +
                                          n +
                                          ",0,1," +
                                          f +
                                          "," +
                                          (t - u) +
                                          "," +
                                          (e - h) +
                                          "A" +
                                          n +
                                          "," +
                                          n +
                                          ",0,1," +
                                          f +
                                          "," +
                                          (this._x1 = d) +
                                          "," +
                                          (this._y1 = g))
                                    : p > a &&
                                      (this._ +=
                                          "A" +
                                          n +
                                          "," +
                                          n +
                                          ",0," +
                                          +(p >= r) +
                                          "," +
                                          f +
                                          "," +
                                          (this._x1 = t + n * Math.cos(l)) +
                                          "," +
                                          (this._y1 = e + n * Math.sin(l)))));
                    },
                    rect: function (t, e, n, r) {
                        this._ +=
                            "M" +
                            (this._x0 = this._x1 = +t) +
                            "," +
                            (this._y0 = this._y1 = +e) +
                            "h" +
                            +n +
                            "v" +
                            +r +
                            "h" +
                            -n +
                            "Z";
                    },
                    toString: function () {
                        return this._;
                    },
                };
                const c = l;
            },
            86: (t, e, n) => {
                "use strict";
                function r(t, e) {
                    return t < e ? -1 : t > e ? 1 : t >= e ? 0 : NaN;
                }
                function i(t) {
                    var e;
                    return (
                        1 === t.length &&
                            ((e = t),
                            (t = function (t, n) {
                                return r(e(t), n);
                            })),
                        {
                            left: function (e, n, r, i) {
                                for (null == r && (r = 0), null == i && (i = e.length); r < i; ) {
                                    var a = (r + i) >>> 1;
                                    t(e[a], n) < 0 ? (r = a + 1) : (i = a);
                                }
                                return r;
                            },
                            right: function (e, n, r, i) {
                                for (null == r && (r = 0), null == i && (i = e.length); r < i; ) {
                                    var a = (r + i) >>> 1;
                                    t(e[a], n) > 0 ? (i = a) : (r = a + 1);
                                }
                                return r;
                            },
                        }
                    );
                }
                (n.r(e),
                    n.d(e, {
                        scaleBand: () => R,
                        scaleDiverging: () => br,
                        scaleDivergingLog: () => xr,
                        scaleDivergingPow: () => $r,
                        scaleDivergingSqrt: () => Ar,
                        scaleDivergingSymlog: () => wr,
                        scaleIdentity: () => mt,
                        scaleImplicit: () => T,
                        scaleLinear: () => pt,
                        scaleLog: () => Ct,
                        scaleOrdinal: () => k,
                        scalePoint: () => L,
                        scalePow: () => Dt,
                        scaleQuantile: () => Et,
                        scaleQuantize: () => Bt,
                        scaleSequential: () => dr,
                        scaleSequentialLog: () => gr,
                        scaleSequentialPow: () => pr,
                        scaleSequentialQuantile: () => vr,
                        scaleSequentialSqrt: () => mr,
                        scaleSequentialSymlog: () => fr,
                        scaleSqrt: () => It,
                        scaleSymlog: () => kt,
                        scaleThreshold: () => Ot,
                        scaleTime: () => lr,
                        scaleUtc: () => cr,
                        tickFormat: () => gt,
                    }));
                var a = i(r),
                    o = a.right;
                a.left;
                const s = o;
                var l = Array.prototype,
                    c = (l.slice, l.map, Math.sqrt(50)),
                    u = Math.sqrt(10),
                    h = Math.sqrt(2);
                function d(t, e, n) {
                    var r,
                        i,
                        a,
                        o,
                        s = -1;
                    if (((n = +n), (t = +t) == (e = +e) && n > 0)) return [t];
                    if (((r = e < t) && ((i = t), (t = e), (e = i)), 0 === (o = g(t, e, n)) || !isFinite(o))) return [];
                    if (o > 0)
                        for (
                            t = Math.ceil(t / o), e = Math.floor(e / o), a = new Array((i = Math.ceil(e - t + 1)));
                            ++s < i;
                        )
                            a[s] = (t + s) * o;
                    else
                        for (
                            t = Math.floor(t * o), e = Math.ceil(e * o), a = new Array((i = Math.ceil(t - e + 1)));
                            ++s < i;
                        )
                            a[s] = (t - s) / o;
                    return (r && a.reverse(), a);
                }
                function g(t, e, n) {
                    var r = (e - t) / Math.max(0, n),
                        i = Math.floor(Math.log(r) / Math.LN10),
                        a = r / Math.pow(10, i);
                    return i >= 0
                        ? (a >= c ? 10 : a >= u ? 5 : a >= h ? 2 : 1) * Math.pow(10, i)
                        : -Math.pow(10, -i) / (a >= c ? 10 : a >= u ? 5 : a >= h ? 2 : 1);
                }
                function f(t, e, n) {
                    var r = Math.abs(e - t) / Math.max(0, n),
                        i = Math.pow(10, Math.floor(Math.log(r) / Math.LN10)),
                        a = r / i;
                    return (a >= c ? (i *= 10) : a >= u ? (i *= 5) : a >= h && (i *= 2), e < t ? -i : i);
                }
                function p(t) {
                    return null === t ? NaN : +t;
                }
                function m(t, e, n) {
                    if ((null == n && (n = p), (r = t.length))) {
                        if ((e = +e) <= 0 || r < 2) return +n(t[0], 0, t);
                        if (e >= 1) return +n(t[r - 1], r - 1, t);
                        var r,
                            i = (r - 1) * e,
                            a = Math.floor(i),
                            o = +n(t[a], a, t);
                        return o + (+n(t[a + 1], a + 1, t) - o) * (i - a);
                    }
                }
                function v(t, e) {
                    switch (arguments.length) {
                        case 0:
                            break;
                        case 1:
                            this.range(t);
                            break;
                        default:
                            this.range(e).domain(t);
                    }
                    return this;
                }
                function y(t, e) {
                    switch (arguments.length) {
                        case 0:
                            break;
                        case 1:
                            this.interpolator(t);
                            break;
                        default:
                            this.interpolator(e).domain(t);
                    }
                    return this;
                }
                var b = "$";
                function x() {}
                function w(t, e) {
                    var n = new x();
                    if (t instanceof x)
                        t.each(function (t, e) {
                            n.set(e, t);
                        });
                    else if (Array.isArray(t)) {
                        var r,
                            i = -1,
                            a = t.length;
                        if (null == e) for (; ++i < a; ) n.set(i, t[i]);
                        else for (; ++i < a; ) n.set(e((r = t[i]), i, t), r);
                    } else if (t) for (var o in t) n.set(o, t[o]);
                    return n;
                }
                x.prototype = w.prototype = {
                    constructor: x,
                    has: function (t) {
                        return b + t in this;
                    },
                    get: function (t) {
                        return this[b + t];
                    },
                    set: function (t, e) {
                        return ((this[b + t] = e), this);
                    },
                    remove: function (t) {
                        var e = b + t;
                        return e in this && delete this[e];
                    },
                    clear: function () {
                        for (var t in this) t[0] === b && delete this[t];
                    },
                    keys: function () {
                        var t = [];
                        for (var e in this) e[0] === b && t.push(e.slice(1));
                        return t;
                    },
                    values: function () {
                        var t = [];
                        for (var e in this) e[0] === b && t.push(this[e]);
                        return t;
                    },
                    entries: function () {
                        var t = [];
                        for (var e in this) e[0] === b && t.push({ key: e.slice(1), value: this[e] });
                        return t;
                    },
                    size: function () {
                        var t = 0;
                        for (var e in this) e[0] === b && ++t;
                        return t;
                    },
                    empty: function () {
                        for (var t in this) if (t[0] === b) return !1;
                        return !0;
                    },
                    each: function (t) {
                        for (var e in this) e[0] === b && t(this[e], e.slice(1), this);
                    },
                };
                const $ = w;
                function A() {}
                var _ = $.prototype;
                A.prototype = function (t, e) {
                    var n = new A();
                    if (t instanceof A)
                        t.each(function (t) {
                            n.add(t);
                        });
                    else if (t) {
                        var r = -1,
                            i = t.length;
                        if (null == e) for (; ++r < i; ) n.add(t[r]);
                        else for (; ++r < i; ) n.add(e(t[r], r, t));
                    }
                    return n;
                }.prototype = {
                    constructor: A,
                    has: _.has,
                    add: function (t) {
                        return ((this[b + (t += "")] = t), this);
                    },
                    remove: _.remove,
                    clear: _.clear,
                    values: _.keys,
                    size: _.size,
                    empty: _.empty,
                    each: _.each,
                };
                var C = Array.prototype,
                    S = C.map,
                    M = C.slice,
                    T = { name: "implicit" };
                function k() {
                    var t = $(),
                        e = [],
                        n = [],
                        r = T;
                    function i(i) {
                        var a = i + "",
                            o = t.get(a);
                        if (!o) {
                            if (r !== T) return r;
                            t.set(a, (o = e.push(i)));
                        }
                        return n[(o - 1) % n.length];
                    }
                    return (
                        (i.domain = function (n) {
                            if (!arguments.length) return e.slice();
                            ((e = []), (t = $()));
                            for (var r, a, o = -1, s = n.length; ++o < s; )
                                t.has((a = (r = n[o]) + "")) || t.set(a, e.push(r));
                            return i;
                        }),
                        (i.range = function (t) {
                            return arguments.length ? ((n = M.call(t)), i) : n.slice();
                        }),
                        (i.unknown = function (t) {
                            return arguments.length ? ((r = t), i) : r;
                        }),
                        (i.copy = function () {
                            return k(e, n).unknown(r);
                        }),
                        v.apply(i, arguments),
                        i
                    );
                }
                function R() {
                    var t,
                        e,
                        n = k().unknown(void 0),
                        r = n.domain,
                        i = n.range,
                        a = [0, 1],
                        o = !1,
                        s = 0,
                        l = 0,
                        c = 0.5;
                    function u() {
                        var n = r().length,
                            u = a[1] < a[0],
                            h = a[u - 0],
                            d = a[1 - u];
                        ((t = (d - h) / Math.max(1, n - s + 2 * l)),
                            o && (t = Math.floor(t)),
                            (h += (d - h - t * (n - s)) * c),
                            (e = t * (1 - s)),
                            o && ((h = Math.round(h)), (e = Math.round(e))));
                        var g = (function (t, e, n) {
                            ((t = +t),
                                (e = +e),
                                (n = (i = arguments.length) < 2 ? ((e = t), (t = 0), 1) : i < 3 ? 1 : +n));
                            for (var r = -1, i = 0 | Math.max(0, Math.ceil((e - t) / n)), a = new Array(i); ++r < i; )
                                a[r] = t + r * n;
                            return a;
                        })(n).map(function (e) {
                            return h + t * e;
                        });
                        return i(u ? g.reverse() : g);
                    }
                    return (
                        delete n.unknown,
                        (n.domain = function (t) {
                            return arguments.length ? (r(t), u()) : r();
                        }),
                        (n.range = function (t) {
                            return arguments.length ? ((a = [+t[0], +t[1]]), u()) : a.slice();
                        }),
                        (n.rangeRound = function (t) {
                            return ((a = [+t[0], +t[1]]), (o = !0), u());
                        }),
                        (n.bandwidth = function () {
                            return e;
                        }),
                        (n.step = function () {
                            return t;
                        }),
                        (n.round = function (t) {
                            return arguments.length ? ((o = !!t), u()) : o;
                        }),
                        (n.padding = function (t) {
                            return arguments.length ? ((s = Math.min(1, (l = +t))), u()) : s;
                        }),
                        (n.paddingInner = function (t) {
                            return arguments.length ? ((s = Math.min(1, t)), u()) : s;
                        }),
                        (n.paddingOuter = function (t) {
                            return arguments.length ? ((l = +t), u()) : l;
                        }),
                        (n.align = function (t) {
                            return arguments.length ? ((c = Math.max(0, Math.min(1, t))), u()) : c;
                        }),
                        (n.copy = function () {
                            return R(r(), a).round(o).paddingInner(s).paddingOuter(l).align(c);
                        }),
                        v.apply(u(), arguments)
                    );
                }
                function N(t) {
                    var e = t.copy;
                    return (
                        (t.padding = t.paddingOuter),
                        delete t.paddingInner,
                        delete t.paddingOuter,
                        (t.copy = function () {
                            return N(e());
                        }),
                        t
                    );
                }
                function L() {
                    return N(R.apply(null, arguments).paddingInner(1));
                }
                var P = n(650),
                    D = n(765);
                function I(t, e) {
                    var n,
                        r = e ? e.length : 0,
                        i = t ? Math.min(r, t.length) : 0,
                        a = new Array(i),
                        o = new Array(r);
                    for (n = 0; n < i; ++n) a[n] = H(t[n], e[n]);
                    for (; n < r; ++n) o[n] = e[n];
                    return function (t) {
                        for (n = 0; n < i; ++n) o[n] = a[n](t);
                        return o;
                    };
                }
                function E(t, e) {
                    var n = new Date();
                    return (
                        (t = +t),
                        (e = +e),
                        function (r) {
                            return (n.setTime(t * (1 - r) + e * r), n);
                        }
                    );
                }
                var B = n(626);
                function O(t, e) {
                    var n,
                        r = {},
                        i = {};
                    for (n in ((null !== t && "object" == typeof t) || (t = {}),
                    (null !== e && "object" == typeof e) || (e = {}),
                    e))
                        n in t ? (r[n] = H(t[n], e[n])) : (i[n] = e[n]);
                    return function (t) {
                        for (n in r) i[n] = r[n](t);
                        return i;
                    };
                }
                var F = n(843),
                    z = n(302);
                function q(t, e) {
                    e || (e = []);
                    var n,
                        r = t ? Math.min(e.length, t.length) : 0,
                        i = e.slice();
                    return function (a) {
                        for (n = 0; n < r; ++n) i[n] = t[n] * (1 - a) + e[n] * a;
                        return i;
                    };
                }
                function H(t, e) {
                    var n,
                        r,
                        i = typeof e;
                    return null == e || "boolean" === i
                        ? (0, z.Z)(e)
                        : ("number" === i
                              ? B.Z
                              : "string" === i
                                ? (n = (0, P.ZP)(e))
                                    ? ((e = n), D.ZP)
                                    : F.Z
                                : e instanceof P.ZP
                                  ? D.ZP
                                  : e instanceof Date
                                    ? E
                                    : ((r = e),
                                      !ArrayBuffer.isView(r) || r instanceof DataView
                                          ? Array.isArray(e)
                                              ? I
                                              : ("function" != typeof e.valueOf && "function" != typeof e.toString) ||
                                                  isNaN(e)
                                                ? O
                                                : B.Z
                                          : q))(t, e);
                }
                function j(t, e) {
                    return (
                        (t = +t),
                        (e = +e),
                        function (n) {
                            return Math.round(t * (1 - n) + e * n);
                        }
                    );
                }
                function U(t) {
                    return +t;
                }
                var V = [0, 1];
                function W(t) {
                    return t;
                }
                function Y(t, e) {
                    return (e -= t = +t)
                        ? function (n) {
                              return (n - t) / e;
                          }
                        : ((n = isNaN(e) ? NaN : 0.5),
                          function () {
                              return n;
                          });
                    var n;
                }
                function Z(t) {
                    var e,
                        n = t[0],
                        r = t[t.length - 1];
                    return (
                        n > r && ((e = n), (n = r), (r = e)),
                        function (t) {
                            return Math.max(n, Math.min(r, t));
                        }
                    );
                }
                function X(t, e, n) {
                    var r = t[0],
                        i = t[1],
                        a = e[0],
                        o = e[1];
                    return (
                        i < r ? ((r = Y(i, r)), (a = n(o, a))) : ((r = Y(r, i)), (a = n(a, o))),
                        function (t) {
                            return a(r(t));
                        }
                    );
                }
                function G(t, e, n) {
                    var r = Math.min(t.length, e.length) - 1,
                        i = new Array(r),
                        a = new Array(r),
                        o = -1;
                    for (t[r] < t[0] && ((t = t.slice().reverse()), (e = e.slice().reverse())); ++o < r; )
                        ((i[o] = Y(t[o], t[o + 1])), (a[o] = n(e[o], e[o + 1])));
                    return function (e) {
                        var n = s(t, e, 1, r) - 1;
                        return a[n](i[n](e));
                    };
                }
                function K(t, e) {
                    return e
                        .domain(t.domain())
                        .range(t.range())
                        .interpolate(t.interpolate())
                        .clamp(t.clamp())
                        .unknown(t.unknown());
                }
                function Q() {
                    var t,
                        e,
                        n,
                        r,
                        i,
                        a,
                        o = V,
                        s = V,
                        l = H,
                        c = W;
                    function u() {
                        return ((r = Math.min(o.length, s.length) > 2 ? G : X), (i = a = null), h);
                    }
                    function h(e) {
                        return isNaN((e = +e)) ? n : (i || (i = r(o.map(t), s, l)))(t(c(e)));
                    }
                    return (
                        (h.invert = function (n) {
                            return c(e((a || (a = r(s, o.map(t), B.Z)))(n)));
                        }),
                        (h.domain = function (t) {
                            return arguments.length ? ((o = S.call(t, U)), c === W || (c = Z(o)), u()) : o.slice();
                        }),
                        (h.range = function (t) {
                            return arguments.length ? ((s = M.call(t)), u()) : s.slice();
                        }),
                        (h.rangeRound = function (t) {
                            return ((s = M.call(t)), (l = j), u());
                        }),
                        (h.clamp = function (t) {
                            return arguments.length ? ((c = t ? Z(o) : W), h) : c !== W;
                        }),
                        (h.interpolate = function (t) {
                            return arguments.length ? ((l = t), u()) : l;
                        }),
                        (h.unknown = function (t) {
                            return arguments.length ? ((n = t), h) : n;
                        }),
                        function (n, r) {
                            return ((t = n), (e = r), u());
                        }
                    );
                }
                function J(t, e) {
                    return Q()(t, e);
                }
                function tt(t, e) {
                    if ((n = (t = e ? t.toExponential(e - 1) : t.toExponential()).indexOf("e")) < 0) return null;
                    var n,
                        r = t.slice(0, n);
                    return [r.length > 1 ? r[0] + r.slice(2) : r, +t.slice(n + 1)];
                }
                function et(t) {
                    return (t = tt(Math.abs(t))) ? t[1] : NaN;
                }
                var nt,
                    rt = /^(?:(.)?([<>=^]))?([+\-( ])?([$#])?(0)?(\d+)?(,)?(\.\d+)?(~)?([a-z%])?$/i;
                function it(t) {
                    return new at(t);
                }
                function at(t) {
                    if (!(e = rt.exec(t))) throw new Error("invalid format: " + t);
                    var e;
                    ((this.fill = e[1] || " "),
                        (this.align = e[2] || ">"),
                        (this.sign = e[3] || "-"),
                        (this.symbol = e[4] || ""),
                        (this.zero = !!e[5]),
                        (this.width = e[6] && +e[6]),
                        (this.comma = !!e[7]),
                        (this.precision = e[8] && +e[8].slice(1)),
                        (this.trim = !!e[9]),
                        (this.type = e[10] || ""));
                }
                function ot(t, e) {
                    var n = tt(t, e);
                    if (!n) return t + "";
                    var r = n[0],
                        i = n[1];
                    return i < 0
                        ? "0." + new Array(-i).join("0") + r
                        : r.length > i + 1
                          ? r.slice(0, i + 1) + "." + r.slice(i + 1)
                          : r + new Array(i - r.length + 2).join("0");
                }
                ((it.prototype = at.prototype),
                    (at.prototype.toString = function () {
                        return (
                            this.fill +
                            this.align +
                            this.sign +
                            this.symbol +
                            (this.zero ? "0" : "") +
                            (null == this.width ? "" : Math.max(1, 0 | this.width)) +
                            (this.comma ? "," : "") +
                            (null == this.precision ? "" : "." + Math.max(0, 0 | this.precision)) +
                            (this.trim ? "~" : "") +
                            this.type
                        );
                    }));
                const st = {
                    "%": function (t, e) {
                        return (100 * t).toFixed(e);
                    },
                    b: function (t) {
                        return Math.round(t).toString(2);
                    },
                    c: function (t) {
                        return t + "";
                    },
                    d: function (t) {
                        return Math.round(t).toString(10);
                    },
                    e: function (t, e) {
                        return t.toExponential(e);
                    },
                    f: function (t, e) {
                        return t.toFixed(e);
                    },
                    g: function (t, e) {
                        return t.toPrecision(e);
                    },
                    o: function (t) {
                        return Math.round(t).toString(8);
                    },
                    p: function (t, e) {
                        return ot(100 * t, e);
                    },
                    r: ot,
                    s: function (t, e) {
                        var n = tt(t, e);
                        if (!n) return t + "";
                        var r = n[0],
                            i = n[1],
                            a = i - (nt = 3 * Math.max(-8, Math.min(8, Math.floor(i / 3)))) + 1,
                            o = r.length;
                        return a === o
                            ? r
                            : a > o
                              ? r + new Array(a - o + 1).join("0")
                              : a > 0
                                ? r.slice(0, a) + "." + r.slice(a)
                                : "0." + new Array(1 - a).join("0") + tt(t, Math.max(0, e + a - 1))[0];
                    },
                    X: function (t) {
                        return Math.round(t).toString(16).toUpperCase();
                    },
                    x: function (t) {
                        return Math.round(t).toString(16);
                    },
                };
                function lt(t) {
                    return t;
                }
                var ct,
                    ut,
                    ht,
                    dt = ["y", "z", "a", "f", "p", "n", "µ", "m", "", "k", "M", "G", "T", "P", "E", "Z", "Y"];
                function gt(t, e, n, r) {
                    var i,
                        a = f(t, e, n);
                    switch ((r = it(null == r ? ",f" : r)).type) {
                        case "s":
                            var o = Math.max(Math.abs(t), Math.abs(e));
                            return (
                                null != r.precision ||
                                    isNaN(
                                        (i = (function (t, e) {
                                            return Math.max(
                                                0,
                                                3 * Math.max(-8, Math.min(8, Math.floor(et(e) / 3))) - et(Math.abs(t)),
                                            );
                                        })(a, o)),
                                    ) ||
                                    (r.precision = i),
                                ht(r, o)
                            );
                        case "":
                        case "e":
                        case "g":
                        case "p":
                        case "r":
                            null != r.precision ||
                                isNaN(
                                    (i = (function (t, e) {
                                        return (
                                            (t = Math.abs(t)),
                                            (e = Math.abs(e) - t),
                                            Math.max(0, et(e) - et(t)) + 1
                                        );
                                    })(a, Math.max(Math.abs(t), Math.abs(e)))),
                                ) ||
                                (r.precision = i - ("e" === r.type));
                            break;
                        case "f":
                        case "%":
                            null != r.precision ||
                                isNaN(
                                    (i = (function (t) {
                                        return Math.max(0, -et(Math.abs(t)));
                                    })(a)),
                                ) ||
                                (r.precision = i - 2 * ("%" === r.type));
                    }
                    return ut(r);
                }
                function ft(t) {
                    var e = t.domain;
                    return (
                        (t.ticks = function (t) {
                            var n = e();
                            return d(n[0], n[n.length - 1], null == t ? 10 : t);
                        }),
                        (t.tickFormat = function (t, n) {
                            var r = e();
                            return gt(r[0], r[r.length - 1], null == t ? 10 : t, n);
                        }),
                        (t.nice = function (n) {
                            null == n && (n = 10);
                            var r,
                                i = e(),
                                a = 0,
                                o = i.length - 1,
                                s = i[a],
                                l = i[o];
                            return (
                                l < s && ((r = s), (s = l), (l = r), (r = a), (a = o), (o = r)),
                                (r = g(s, l, n)) > 0
                                    ? (r = g((s = Math.floor(s / r) * r), (l = Math.ceil(l / r) * r), n))
                                    : r < 0 && (r = g((s = Math.ceil(s * r) / r), (l = Math.floor(l * r) / r), n)),
                                r > 0
                                    ? ((i[a] = Math.floor(s / r) * r), (i[o] = Math.ceil(l / r) * r), e(i))
                                    : r < 0 && ((i[a] = Math.ceil(s * r) / r), (i[o] = Math.floor(l * r) / r), e(i)),
                                t
                            );
                        }),
                        t
                    );
                }
                function pt() {
                    var t = J(W, W);
                    return (
                        (t.copy = function () {
                            return K(t, pt());
                        }),
                        v.apply(t, arguments),
                        ft(t)
                    );
                }
                function mt(t) {
                    var e;
                    function n(t) {
                        return isNaN((t = +t)) ? e : t;
                    }
                    return (
                        (n.invert = n),
                        (n.domain = n.range =
                            function (e) {
                                return arguments.length ? ((t = S.call(e, U)), n) : t.slice();
                            }),
                        (n.unknown = function (t) {
                            return arguments.length ? ((e = t), n) : e;
                        }),
                        (n.copy = function () {
                            return mt(t).unknown(e);
                        }),
                        (t = arguments.length ? S.call(t, U) : [0, 1]),
                        ft(n)
                    );
                }
                function vt(t, e) {
                    var n,
                        r = 0,
                        i = (t = t.slice()).length - 1,
                        a = t[r],
                        o = t[i];
                    return (
                        o < a && ((n = r), (r = i), (i = n), (n = a), (a = o), (o = n)),
                        (t[r] = e.floor(a)),
                        (t[i] = e.ceil(o)),
                        t
                    );
                }
                function yt(t) {
                    return Math.log(t);
                }
                function bt(t) {
                    return Math.exp(t);
                }
                function xt(t) {
                    return -Math.log(-t);
                }
                function wt(t) {
                    return -Math.exp(-t);
                }
                function $t(t) {
                    return isFinite(t) ? +("1e" + t) : t < 0 ? 0 : t;
                }
                function At(t) {
                    return function (e) {
                        return -t(-e);
                    };
                }
                function _t(t) {
                    var e,
                        n,
                        r = t(yt, bt),
                        i = r.domain,
                        a = 10;
                    function o() {
                        return (
                            (e = (function (t) {
                                return t === Math.E
                                    ? Math.log
                                    : (10 === t && Math.log10) ||
                                          (2 === t && Math.log2) ||
                                          ((t = Math.log(t)),
                                          function (e) {
                                              return Math.log(e) / t;
                                          });
                            })(a)),
                            (n = (function (t) {
                                return 10 === t
                                    ? $t
                                    : t === Math.E
                                      ? Math.exp
                                      : function (e) {
                                            return Math.pow(t, e);
                                        };
                            })(a)),
                            i()[0] < 0 ? ((e = At(e)), (n = At(n)), t(xt, wt)) : t(yt, bt),
                            r
                        );
                    }
                    return (
                        (r.base = function (t) {
                            return arguments.length ? ((a = +t), o()) : a;
                        }),
                        (r.domain = function (t) {
                            return arguments.length ? (i(t), o()) : i();
                        }),
                        (r.ticks = function (t) {
                            var r,
                                o = i(),
                                s = o[0],
                                l = o[o.length - 1];
                            (r = l < s) && ((g = s), (s = l), (l = g));
                            var c,
                                u,
                                h,
                                g = e(s),
                                f = e(l),
                                p = null == t ? 10 : +t,
                                m = [];
                            if (!(a % 1) && f - g < p) {
                                if (((g = Math.round(g) - 1), (f = Math.round(f) + 1), s > 0)) {
                                    for (; g < f; ++g)
                                        for (u = 1, c = n(g); u < a; ++u)
                                            if (!((h = c * u) < s)) {
                                                if (h > l) break;
                                                m.push(h);
                                            }
                                } else
                                    for (; g < f; ++g)
                                        for (u = a - 1, c = n(g); u >= 1; --u)
                                            if (!((h = c * u) < s)) {
                                                if (h > l) break;
                                                m.push(h);
                                            }
                            } else m = d(g, f, Math.min(f - g, p)).map(n);
                            return r ? m.reverse() : m;
                        }),
                        (r.tickFormat = function (t, i) {
                            if (
                                (null == i && (i = 10 === a ? ".0e" : ","),
                                "function" != typeof i && (i = ut(i)),
                                t === 1 / 0)
                            )
                                return i;
                            null == t && (t = 10);
                            var o = Math.max(1, (a * t) / r.ticks().length);
                            return function (t) {
                                var r = t / n(Math.round(e(t)));
                                return (r * a < a - 0.5 && (r *= a), r <= o ? i(t) : "");
                            };
                        }),
                        (r.nice = function () {
                            return i(
                                vt(i(), {
                                    floor: function (t) {
                                        return n(Math.floor(e(t)));
                                    },
                                    ceil: function (t) {
                                        return n(Math.ceil(e(t)));
                                    },
                                }),
                            );
                        }),
                        r
                    );
                }
                function Ct() {
                    var t = _t(Q()).domain([1, 10]);
                    return (
                        (t.copy = function () {
                            return K(t, Ct()).base(t.base());
                        }),
                        v.apply(t, arguments),
                        t
                    );
                }
                function St(t) {
                    return function (e) {
                        return Math.sign(e) * Math.log1p(Math.abs(e / t));
                    };
                }
                function Mt(t) {
                    return function (e) {
                        return Math.sign(e) * Math.expm1(Math.abs(e)) * t;
                    };
                }
                function Tt(t) {
                    var e = 1,
                        n = t(St(e), Mt(e));
                    return (
                        (n.constant = function (n) {
                            return arguments.length ? t(St((e = +n)), Mt(e)) : e;
                        }),
                        ft(n)
                    );
                }
                function kt() {
                    var t = Tt(Q());
                    return (
                        (t.copy = function () {
                            return K(t, kt()).constant(t.constant());
                        }),
                        v.apply(t, arguments)
                    );
                }
                function Rt(t) {
                    return function (e) {
                        return e < 0 ? -Math.pow(-e, t) : Math.pow(e, t);
                    };
                }
                function Nt(t) {
                    return t < 0 ? -Math.sqrt(-t) : Math.sqrt(t);
                }
                function Lt(t) {
                    return t < 0 ? -t * t : t * t;
                }
                function Pt(t) {
                    var e = t(W, W),
                        n = 1;
                    return (
                        (e.exponent = function (e) {
                            return arguments.length
                                ? 1 == (n = +e)
                                    ? t(W, W)
                                    : 0.5 === n
                                      ? t(Nt, Lt)
                                      : t(Rt(n), Rt(1 / n))
                                : n;
                        }),
                        ft(e)
                    );
                }
                function Dt() {
                    var t = Pt(Q());
                    return (
                        (t.copy = function () {
                            return K(t, Dt()).exponent(t.exponent());
                        }),
                        v.apply(t, arguments),
                        t
                    );
                }
                function It() {
                    return Dt.apply(null, arguments).exponent(0.5);
                }
                function Et() {
                    var t,
                        e = [],
                        n = [],
                        i = [];
                    function a() {
                        var t = 0,
                            r = Math.max(1, n.length);
                        for (i = new Array(r - 1); ++t < r; ) i[t - 1] = m(e, t / r);
                        return o;
                    }
                    function o(e) {
                        return isNaN((e = +e)) ? t : n[s(i, e)];
                    }
                    return (
                        (o.invertExtent = function (t) {
                            var r = n.indexOf(t);
                            return r < 0
                                ? [NaN, NaN]
                                : [r > 0 ? i[r - 1] : e[0], r < i.length ? i[r] : e[e.length - 1]];
                        }),
                        (o.domain = function (t) {
                            if (!arguments.length) return e.slice();
                            e = [];
                            for (var n, i = 0, o = t.length; i < o; ++i)
                                null == (n = t[i]) || isNaN((n = +n)) || e.push(n);
                            return (e.sort(r), a());
                        }),
                        (o.range = function (t) {
                            return arguments.length ? ((n = M.call(t)), a()) : n.slice();
                        }),
                        (o.unknown = function (e) {
                            return arguments.length ? ((t = e), o) : t;
                        }),
                        (o.quantiles = function () {
                            return i.slice();
                        }),
                        (o.copy = function () {
                            return Et().domain(e).range(n).unknown(t);
                        }),
                        v.apply(o, arguments)
                    );
                }
                function Bt() {
                    var t,
                        e = 0,
                        n = 1,
                        r = 1,
                        i = [0.5],
                        a = [0, 1];
                    function o(e) {
                        return e <= e ? a[s(i, e, 0, r)] : t;
                    }
                    function l() {
                        var t = -1;
                        for (i = new Array(r); ++t < r; ) i[t] = ((t + 1) * n - (t - r) * e) / (r + 1);
                        return o;
                    }
                    return (
                        (o.domain = function (t) {
                            return arguments.length ? ((e = +t[0]), (n = +t[1]), l()) : [e, n];
                        }),
                        (o.range = function (t) {
                            return arguments.length ? ((r = (a = M.call(t)).length - 1), l()) : a.slice();
                        }),
                        (o.invertExtent = function (t) {
                            var o = a.indexOf(t);
                            return o < 0 ? [NaN, NaN] : o < 1 ? [e, i[0]] : o >= r ? [i[r - 1], n] : [i[o - 1], i[o]];
                        }),
                        (o.unknown = function (e) {
                            return arguments.length ? ((t = e), o) : o;
                        }),
                        (o.thresholds = function () {
                            return i.slice();
                        }),
                        (o.copy = function () {
                            return Bt().domain([e, n]).range(a).unknown(t);
                        }),
                        v.apply(ft(o), arguments)
                    );
                }
                function Ot() {
                    var t,
                        e = [0.5],
                        n = [0, 1],
                        r = 1;
                    function i(i) {
                        return i <= i ? n[s(e, i, 0, r)] : t;
                    }
                    return (
                        (i.domain = function (t) {
                            return arguments.length
                                ? ((e = M.call(t)), (r = Math.min(e.length, n.length - 1)), i)
                                : e.slice();
                        }),
                        (i.range = function (t) {
                            return arguments.length
                                ? ((n = M.call(t)), (r = Math.min(e.length, n.length - 1)), i)
                                : n.slice();
                        }),
                        (i.invertExtent = function (t) {
                            var r = n.indexOf(t);
                            return [e[r - 1], e[r]];
                        }),
                        (i.unknown = function (e) {
                            return arguments.length ? ((t = e), i) : t;
                        }),
                        (i.copy = function () {
                            return Ot().domain(e).range(n).unknown(t);
                        }),
                        v.apply(i, arguments)
                    );
                }
                ((ct = (function (t) {
                    var e,
                        n,
                        r =
                            t.grouping && t.thousands
                                ? ((e = t.grouping),
                                  (n = t.thousands),
                                  function (t, r) {
                                      for (
                                          var i = t.length, a = [], o = 0, s = e[0], l = 0;
                                          i > 0 &&
                                          s > 0 &&
                                          (l + s + 1 > r && (s = Math.max(1, r - l)),
                                          a.push(t.substring((i -= s), i + s)),
                                          !((l += s + 1) > r));
                                      )
                                          s = e[(o = (o + 1) % e.length)];
                                      return a.reverse().join(n);
                                  })
                                : lt,
                        i = t.currency,
                        a = t.decimal,
                        o = t.numerals
                            ? (function (t) {
                                  return function (e) {
                                      return e.replace(/[0-9]/g, function (e) {
                                          return t[+e];
                                      });
                                  };
                              })(t.numerals)
                            : lt,
                        s = t.percent || "%";
                    function l(t) {
                        var e = (t = it(t)).fill,
                            n = t.align,
                            l = t.sign,
                            c = t.symbol,
                            u = t.zero,
                            h = t.width,
                            d = t.comma,
                            g = t.precision,
                            f = t.trim,
                            p = t.type;
                        ("n" === p ? ((d = !0), (p = "g")) : st[p] || (null == g && (g = 12), (f = !0), (p = "g")),
                            (u || ("0" === e && "=" === n)) && ((u = !0), (e = "0"), (n = "=")));
                        var m = "$" === c ? i[0] : "#" === c && /[boxX]/.test(p) ? "0" + p.toLowerCase() : "",
                            v = "$" === c ? i[1] : /[%p]/.test(p) ? s : "",
                            y = st[p],
                            b = /[defgprs%]/.test(p);
                        function x(t) {
                            var i,
                                s,
                                c,
                                x = m,
                                w = v;
                            if ("c" === p) ((w = y(t) + w), (t = ""));
                            else {
                                var $ = (t = +t) < 0;
                                if (
                                    ((t = y(Math.abs(t), g)),
                                    f &&
                                        (t = (function (t) {
                                            t: for (var e, n = t.length, r = 1, i = -1; r < n; ++r)
                                                switch (t[r]) {
                                                    case ".":
                                                        i = e = r;
                                                        break;
                                                    case "0":
                                                        (0 === i && (i = r), (e = r));
                                                        break;
                                                    default:
                                                        if (i > 0) {
                                                            if (!+t[r]) break t;
                                                            i = 0;
                                                        }
                                                }
                                            return i > 0 ? t.slice(0, i) + t.slice(e + 1) : t;
                                        })(t)),
                                    $ && 0 == +t && ($ = !1),
                                    (x = ($ ? ("(" === l ? l : "-") : "-" === l || "(" === l ? "" : l) + x),
                                    (w = ("s" === p ? dt[8 + nt / 3] : "") + w + ($ && "(" === l ? ")" : "")),
                                    b)
                                )
                                    for (i = -1, s = t.length; ++i < s; )
                                        if (48 > (c = t.charCodeAt(i)) || c > 57) {
                                            ((w = (46 === c ? a + t.slice(i + 1) : t.slice(i)) + w),
                                                (t = t.slice(0, i)));
                                            break;
                                        }
                            }
                            d && !u && (t = r(t, 1 / 0));
                            var A = x.length + t.length + w.length,
                                _ = A < h ? new Array(h - A + 1).join(e) : "";
                            switch ((d && u && ((t = r(_ + t, _.length ? h - w.length : 1 / 0)), (_ = "")), n)) {
                                case "<":
                                    t = x + t + w + _;
                                    break;
                                case "=":
                                    t = x + _ + t + w;
                                    break;
                                case "^":
                                    t = _.slice(0, (A = _.length >> 1)) + x + t + w + _.slice(A);
                                    break;
                                default:
                                    t = _ + x + t + w;
                            }
                            return o(t);
                        }
                        return (
                            (g =
                                null == g
                                    ? 6
                                    : /[gprs]/.test(p)
                                      ? Math.max(1, Math.min(21, g))
                                      : Math.max(0, Math.min(20, g))),
                            (x.toString = function () {
                                return t + "";
                            }),
                            x
                        );
                    }
                    return {
                        format: l,
                        formatPrefix: function (t, e) {
                            var n = l((((t = it(t)).type = "f"), t)),
                                r = 3 * Math.max(-8, Math.min(8, Math.floor(et(e) / 3))),
                                i = Math.pow(10, -r),
                                a = dt[8 + r / 3];
                            return function (t) {
                                return n(i * t) + a;
                            };
                        },
                    };
                })({ decimal: ".", thousands: ",", grouping: [3], currency: ["$", ""] })),
                    (ut = ct.format),
                    (ht = ct.formatPrefix));
                var Ft = new Date(),
                    zt = new Date();
                function qt(t, e, n, r) {
                    function i(e) {
                        return (t((e = new Date(+e))), e);
                    }
                    return (
                        (i.floor = i),
                        (i.ceil = function (n) {
                            return (t((n = new Date(n - 1))), e(n, 1), t(n), n);
                        }),
                        (i.round = function (t) {
                            var e = i(t),
                                n = i.ceil(t);
                            return t - e < n - t ? e : n;
                        }),
                        (i.offset = function (t, n) {
                            return (e((t = new Date(+t)), null == n ? 1 : Math.floor(n)), t);
                        }),
                        (i.range = function (n, r, a) {
                            var o,
                                s = [];
                            if (((n = i.ceil(n)), (a = null == a ? 1 : Math.floor(a)), !(n < r && a > 0))) return s;
                            do {
                                (s.push((o = new Date(+n))), e(n, a), t(n));
                            } while (o < n && n < r);
                            return s;
                        }),
                        (i.filter = function (n) {
                            return qt(
                                function (e) {
                                    if (e >= e) for (; t(e), !n(e); ) e.setTime(e - 1);
                                },
                                function (t, r) {
                                    if (t >= t)
                                        if (r < 0) for (; ++r <= 0; ) for (; e(t, -1), !n(t); );
                                        else for (; --r >= 0; ) for (; e(t, 1), !n(t); );
                                },
                            );
                        }),
                        n &&
                            ((i.count = function (e, r) {
                                return (Ft.setTime(+e), zt.setTime(+r), t(Ft), t(zt), Math.floor(n(Ft, zt)));
                            }),
                            (i.every = function (t) {
                                return (
                                    (t = Math.floor(t)),
                                    isFinite(t) && t > 0
                                        ? t > 1
                                            ? i.filter(
                                                  r
                                                      ? function (e) {
                                                            return r(e) % t == 0;
                                                        }
                                                      : function (e) {
                                                            return i.count(0, e) % t == 0;
                                                        },
                                              )
                                            : i
                                        : null
                                );
                            })),
                        i
                    );
                }
                var Ht = qt(
                    function () {},
                    function (t, e) {
                        t.setTime(+t + e);
                    },
                    function (t, e) {
                        return e - t;
                    },
                );
                Ht.every = function (t) {
                    return (
                        (t = Math.floor(t)),
                        isFinite(t) && t > 0
                            ? t > 1
                                ? qt(
                                      function (e) {
                                          e.setTime(Math.floor(e / t) * t);
                                      },
                                      function (e, n) {
                                          e.setTime(+e + n * t);
                                      },
                                      function (e, n) {
                                          return (n - e) / t;
                                      },
                                  )
                                : Ht
                            : null
                    );
                };
                const jt = Ht;
                Ht.range;
                var Ut = 1e3,
                    Vt = 6e4,
                    Wt = 36e5,
                    Yt = 864e5,
                    Zt = 6048e5,
                    Xt = qt(
                        function (t) {
                            t.setTime(Math.floor(t / Ut) * Ut);
                        },
                        function (t, e) {
                            t.setTime(+t + e * Ut);
                        },
                        function (t, e) {
                            return (e - t) / Ut;
                        },
                        function (t) {
                            return t.getUTCSeconds();
                        },
                    );
                const Gt = Xt;
                Xt.range;
                var Kt = qt(
                    function (t) {
                        t.setTime(Math.floor(t / Vt) * Vt);
                    },
                    function (t, e) {
                        t.setTime(+t + e * Vt);
                    },
                    function (t, e) {
                        return (e - t) / Vt;
                    },
                    function (t) {
                        return t.getMinutes();
                    },
                );
                const Qt = Kt;
                Kt.range;
                var Jt = qt(
                    function (t) {
                        var e = (t.getTimezoneOffset() * Vt) % Wt;
                        (e < 0 && (e += Wt), t.setTime(Math.floor((+t - e) / Wt) * Wt + e));
                    },
                    function (t, e) {
                        t.setTime(+t + e * Wt);
                    },
                    function (t, e) {
                        return (e - t) / Wt;
                    },
                    function (t) {
                        return t.getHours();
                    },
                );
                const te = Jt;
                Jt.range;
                var ee = qt(
                    function (t) {
                        t.setHours(0, 0, 0, 0);
                    },
                    function (t, e) {
                        t.setDate(t.getDate() + e);
                    },
                    function (t, e) {
                        return (e - t - (e.getTimezoneOffset() - t.getTimezoneOffset()) * Vt) / Yt;
                    },
                    function (t) {
                        return t.getDate() - 1;
                    },
                );
                const ne = ee;
                function re(t) {
                    return qt(
                        function (e) {
                            (e.setDate(e.getDate() - ((e.getDay() + 7 - t) % 7)), e.setHours(0, 0, 0, 0));
                        },
                        function (t, e) {
                            t.setDate(t.getDate() + 7 * e);
                        },
                        function (t, e) {
                            return (e - t - (e.getTimezoneOffset() - t.getTimezoneOffset()) * Vt) / Zt;
                        },
                    );
                }
                ee.range;
                var ie = re(0),
                    ae = re(1),
                    oe = re(2),
                    se = re(3),
                    le = re(4),
                    ce = re(5),
                    ue = re(6),
                    he =
                        (ie.range,
                        ae.range,
                        oe.range,
                        se.range,
                        le.range,
                        ce.range,
                        ue.range,
                        qt(
                            function (t) {
                                (t.setDate(1), t.setHours(0, 0, 0, 0));
                            },
                            function (t, e) {
                                t.setMonth(t.getMonth() + e);
                            },
                            function (t, e) {
                                return e.getMonth() - t.getMonth() + 12 * (e.getFullYear() - t.getFullYear());
                            },
                            function (t) {
                                return t.getMonth();
                            },
                        ));
                const de = he;
                he.range;
                var ge = qt(
                    function (t) {
                        (t.setMonth(0, 1), t.setHours(0, 0, 0, 0));
                    },
                    function (t, e) {
                        t.setFullYear(t.getFullYear() + e);
                    },
                    function (t, e) {
                        return e.getFullYear() - t.getFullYear();
                    },
                    function (t) {
                        return t.getFullYear();
                    },
                );
                ge.every = function (t) {
                    return isFinite((t = Math.floor(t))) && t > 0
                        ? qt(
                              function (e) {
                                  (e.setFullYear(Math.floor(e.getFullYear() / t) * t),
                                      e.setMonth(0, 1),
                                      e.setHours(0, 0, 0, 0));
                              },
                              function (e, n) {
                                  e.setFullYear(e.getFullYear() + n * t);
                              },
                          )
                        : null;
                };
                const fe = ge;
                ge.range;
                var pe = qt(
                    function (t) {
                        t.setUTCSeconds(0, 0);
                    },
                    function (t, e) {
                        t.setTime(+t + e * Vt);
                    },
                    function (t, e) {
                        return (e - t) / Vt;
                    },
                    function (t) {
                        return t.getUTCMinutes();
                    },
                );
                const me = pe;
                pe.range;
                var ve = qt(
                    function (t) {
                        t.setUTCMinutes(0, 0, 0);
                    },
                    function (t, e) {
                        t.setTime(+t + e * Wt);
                    },
                    function (t, e) {
                        return (e - t) / Wt;
                    },
                    function (t) {
                        return t.getUTCHours();
                    },
                );
                const ye = ve;
                ve.range;
                var be = qt(
                    function (t) {
                        t.setUTCHours(0, 0, 0, 0);
                    },
                    function (t, e) {
                        t.setUTCDate(t.getUTCDate() + e);
                    },
                    function (t, e) {
                        return (e - t) / Yt;
                    },
                    function (t) {
                        return t.getUTCDate() - 1;
                    },
                );
                const xe = be;
                function we(t) {
                    return qt(
                        function (e) {
                            (e.setUTCDate(e.getUTCDate() - ((e.getUTCDay() + 7 - t) % 7)), e.setUTCHours(0, 0, 0, 0));
                        },
                        function (t, e) {
                            t.setUTCDate(t.getUTCDate() + 7 * e);
                        },
                        function (t, e) {
                            return (e - t) / Zt;
                        },
                    );
                }
                be.range;
                var $e = we(0),
                    Ae = we(1),
                    _e = we(2),
                    Ce = we(3),
                    Se = we(4),
                    Me = we(5),
                    Te = we(6),
                    ke =
                        ($e.range,
                        Ae.range,
                        _e.range,
                        Ce.range,
                        Se.range,
                        Me.range,
                        Te.range,
                        qt(
                            function (t) {
                                (t.setUTCDate(1), t.setUTCHours(0, 0, 0, 0));
                            },
                            function (t, e) {
                                t.setUTCMonth(t.getUTCMonth() + e);
                            },
                            function (t, e) {
                                return (
                                    e.getUTCMonth() - t.getUTCMonth() + 12 * (e.getUTCFullYear() - t.getUTCFullYear())
                                );
                            },
                            function (t) {
                                return t.getUTCMonth();
                            },
                        ));
                const Re = ke;
                ke.range;
                var Ne = qt(
                    function (t) {
                        (t.setUTCMonth(0, 1), t.setUTCHours(0, 0, 0, 0));
                    },
                    function (t, e) {
                        t.setUTCFullYear(t.getUTCFullYear() + e);
                    },
                    function (t, e) {
                        return e.getUTCFullYear() - t.getUTCFullYear();
                    },
                    function (t) {
                        return t.getUTCFullYear();
                    },
                );
                Ne.every = function (t) {
                    return isFinite((t = Math.floor(t))) && t > 0
                        ? qt(
                              function (e) {
                                  (e.setUTCFullYear(Math.floor(e.getUTCFullYear() / t) * t),
                                      e.setUTCMonth(0, 1),
                                      e.setUTCHours(0, 0, 0, 0));
                              },
                              function (e, n) {
                                  e.setUTCFullYear(e.getUTCFullYear() + n * t);
                              },
                          )
                        : null;
                };
                const Le = Ne;
                function Pe(t) {
                    if (0 <= t.y && t.y < 100) {
                        var e = new Date(-1, t.m, t.d, t.H, t.M, t.S, t.L);
                        return (e.setFullYear(t.y), e);
                    }
                    return new Date(t.y, t.m, t.d, t.H, t.M, t.S, t.L);
                }
                function De(t) {
                    if (0 <= t.y && t.y < 100) {
                        var e = new Date(Date.UTC(-1, t.m, t.d, t.H, t.M, t.S, t.L));
                        return (e.setUTCFullYear(t.y), e);
                    }
                    return new Date(Date.UTC(t.y, t.m, t.d, t.H, t.M, t.S, t.L));
                }
                function Ie(t) {
                    return { y: t, m: 0, d: 1, H: 0, M: 0, S: 0, L: 0 };
                }
                Ne.range;
                var Ee,
                    Be,
                    Oe,
                    Fe,
                    ze = { "-": "", _: " ", 0: "0" },
                    qe = /^\s*\d+/,
                    He = /^%/,
                    je = /[\\^$*+?|[\]().{}]/g;
                function Ue(t, e, n) {
                    var r = t < 0 ? "-" : "",
                        i = (r ? -t : t) + "",
                        a = i.length;
                    return r + (a < n ? new Array(n - a + 1).join(e) + i : i);
                }
                function Ve(t) {
                    return t.replace(je, "\\$&");
                }
                function We(t) {
                    return new RegExp("^(?:" + t.map(Ve).join("|") + ")", "i");
                }
                function Ye(t) {
                    for (var e = {}, n = -1, r = t.length; ++n < r; ) e[t[n].toLowerCase()] = n;
                    return e;
                }
                function Ze(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 1));
                    return r ? ((t.w = +r[0]), n + r[0].length) : -1;
                }
                function Xe(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 1));
                    return r ? ((t.u = +r[0]), n + r[0].length) : -1;
                }
                function Ge(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.U = +r[0]), n + r[0].length) : -1;
                }
                function Ke(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.V = +r[0]), n + r[0].length) : -1;
                }
                function Qe(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.W = +r[0]), n + r[0].length) : -1;
                }
                function Je(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 4));
                    return r ? ((t.y = +r[0]), n + r[0].length) : -1;
                }
                function tn(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.y = +r[0] + (+r[0] > 68 ? 1900 : 2e3)), n + r[0].length) : -1;
                }
                function en(t, e, n) {
                    var r = /^(Z)|([+-]\d\d)(?::?(\d\d))?/.exec(e.slice(n, n + 6));
                    return r ? ((t.Z = r[1] ? 0 : -(r[2] + (r[3] || "00"))), n + r[0].length) : -1;
                }
                function nn(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.m = r[0] - 1), n + r[0].length) : -1;
                }
                function rn(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.d = +r[0]), n + r[0].length) : -1;
                }
                function an(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 3));
                    return r ? ((t.m = 0), (t.d = +r[0]), n + r[0].length) : -1;
                }
                function on(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.H = +r[0]), n + r[0].length) : -1;
                }
                function sn(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.M = +r[0]), n + r[0].length) : -1;
                }
                function ln(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 2));
                    return r ? ((t.S = +r[0]), n + r[0].length) : -1;
                }
                function cn(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 3));
                    return r ? ((t.L = +r[0]), n + r[0].length) : -1;
                }
                function un(t, e, n) {
                    var r = qe.exec(e.slice(n, n + 6));
                    return r ? ((t.L = Math.floor(r[0] / 1e3)), n + r[0].length) : -1;
                }
                function hn(t, e, n) {
                    var r = He.exec(e.slice(n, n + 1));
                    return r ? n + r[0].length : -1;
                }
                function dn(t, e, n) {
                    var r = qe.exec(e.slice(n));
                    return r ? ((t.Q = +r[0]), n + r[0].length) : -1;
                }
                function gn(t, e, n) {
                    var r = qe.exec(e.slice(n));
                    return r ? ((t.Q = 1e3 * +r[0]), n + r[0].length) : -1;
                }
                function fn(t, e) {
                    return Ue(t.getDate(), e, 2);
                }
                function pn(t, e) {
                    return Ue(t.getHours(), e, 2);
                }
                function mn(t, e) {
                    return Ue(t.getHours() % 12 || 12, e, 2);
                }
                function vn(t, e) {
                    return Ue(1 + ne.count(fe(t), t), e, 3);
                }
                function yn(t, e) {
                    return Ue(t.getMilliseconds(), e, 3);
                }
                function bn(t, e) {
                    return yn(t, e) + "000";
                }
                function xn(t, e) {
                    return Ue(t.getMonth() + 1, e, 2);
                }
                function wn(t, e) {
                    return Ue(t.getMinutes(), e, 2);
                }
                function $n(t, e) {
                    return Ue(t.getSeconds(), e, 2);
                }
                function An(t) {
                    var e = t.getDay();
                    return 0 === e ? 7 : e;
                }
                function _n(t, e) {
                    return Ue(ie.count(fe(t), t), e, 2);
                }
                function Cn(t, e) {
                    var n = t.getDay();
                    return (
                        (t = n >= 4 || 0 === n ? le(t) : le.ceil(t)),
                        Ue(le.count(fe(t), t) + (4 === fe(t).getDay()), e, 2)
                    );
                }
                function Sn(t) {
                    return t.getDay();
                }
                function Mn(t, e) {
                    return Ue(ae.count(fe(t), t), e, 2);
                }
                function Tn(t, e) {
                    return Ue(t.getFullYear() % 100, e, 2);
                }
                function kn(t, e) {
                    return Ue(t.getFullYear() % 1e4, e, 4);
                }
                function Rn(t) {
                    var e = t.getTimezoneOffset();
                    return (e > 0 ? "-" : ((e *= -1), "+")) + Ue((e / 60) | 0, "0", 2) + Ue(e % 60, "0", 2);
                }
                function Nn(t, e) {
                    return Ue(t.getUTCDate(), e, 2);
                }
                function Ln(t, e) {
                    return Ue(t.getUTCHours(), e, 2);
                }
                function Pn(t, e) {
                    return Ue(t.getUTCHours() % 12 || 12, e, 2);
                }
                function Dn(t, e) {
                    return Ue(1 + xe.count(Le(t), t), e, 3);
                }
                function In(t, e) {
                    return Ue(t.getUTCMilliseconds(), e, 3);
                }
                function En(t, e) {
                    return In(t, e) + "000";
                }
                function Bn(t, e) {
                    return Ue(t.getUTCMonth() + 1, e, 2);
                }
                function On(t, e) {
                    return Ue(t.getUTCMinutes(), e, 2);
                }
                function Fn(t, e) {
                    return Ue(t.getUTCSeconds(), e, 2);
                }
                function zn(t) {
                    var e = t.getUTCDay();
                    return 0 === e ? 7 : e;
                }
                function qn(t, e) {
                    return Ue($e.count(Le(t), t), e, 2);
                }
                function Hn(t, e) {
                    var n = t.getUTCDay();
                    return (
                        (t = n >= 4 || 0 === n ? Se(t) : Se.ceil(t)),
                        Ue(Se.count(Le(t), t) + (4 === Le(t).getUTCDay()), e, 2)
                    );
                }
                function jn(t) {
                    return t.getUTCDay();
                }
                function Un(t, e) {
                    return Ue(Ae.count(Le(t), t), e, 2);
                }
                function Vn(t, e) {
                    return Ue(t.getUTCFullYear() % 100, e, 2);
                }
                function Wn(t, e) {
                    return Ue(t.getUTCFullYear() % 1e4, e, 4);
                }
                function Yn() {
                    return "+0000";
                }
                function Zn() {
                    return "%";
                }
                function Xn(t) {
                    return +t;
                }
                function Gn(t) {
                    return Math.floor(+t / 1e3);
                }
                ((Ee = (function (t) {
                    var e = t.dateTime,
                        n = t.date,
                        r = t.time,
                        i = t.periods,
                        a = t.days,
                        o = t.shortDays,
                        s = t.months,
                        l = t.shortMonths,
                        c = We(i),
                        u = Ye(i),
                        h = We(a),
                        d = Ye(a),
                        g = We(o),
                        f = Ye(o),
                        p = We(s),
                        m = Ye(s),
                        v = We(l),
                        y = Ye(l),
                        b = {
                            a: function (t) {
                                return o[t.getDay()];
                            },
                            A: function (t) {
                                return a[t.getDay()];
                            },
                            b: function (t) {
                                return l[t.getMonth()];
                            },
                            B: function (t) {
                                return s[t.getMonth()];
                            },
                            c: null,
                            d: fn,
                            e: fn,
                            f: bn,
                            H: pn,
                            I: mn,
                            j: vn,
                            L: yn,
                            m: xn,
                            M: wn,
                            p: function (t) {
                                return i[+(t.getHours() >= 12)];
                            },
                            Q: Xn,
                            s: Gn,
                            S: $n,
                            u: An,
                            U: _n,
                            V: Cn,
                            w: Sn,
                            W: Mn,
                            x: null,
                            X: null,
                            y: Tn,
                            Y: kn,
                            Z: Rn,
                            "%": Zn,
                        },
                        x = {
                            a: function (t) {
                                return o[t.getUTCDay()];
                            },
                            A: function (t) {
                                return a[t.getUTCDay()];
                            },
                            b: function (t) {
                                return l[t.getUTCMonth()];
                            },
                            B: function (t) {
                                return s[t.getUTCMonth()];
                            },
                            c: null,
                            d: Nn,
                            e: Nn,
                            f: En,
                            H: Ln,
                            I: Pn,
                            j: Dn,
                            L: In,
                            m: Bn,
                            M: On,
                            p: function (t) {
                                return i[+(t.getUTCHours() >= 12)];
                            },
                            Q: Xn,
                            s: Gn,
                            S: Fn,
                            u: zn,
                            U: qn,
                            V: Hn,
                            w: jn,
                            W: Un,
                            x: null,
                            X: null,
                            y: Vn,
                            Y: Wn,
                            Z: Yn,
                            "%": Zn,
                        },
                        w = {
                            a: function (t, e, n) {
                                var r = g.exec(e.slice(n));
                                return r ? ((t.w = f[r[0].toLowerCase()]), n + r[0].length) : -1;
                            },
                            A: function (t, e, n) {
                                var r = h.exec(e.slice(n));
                                return r ? ((t.w = d[r[0].toLowerCase()]), n + r[0].length) : -1;
                            },
                            b: function (t, e, n) {
                                var r = v.exec(e.slice(n));
                                return r ? ((t.m = y[r[0].toLowerCase()]), n + r[0].length) : -1;
                            },
                            B: function (t, e, n) {
                                var r = p.exec(e.slice(n));
                                return r ? ((t.m = m[r[0].toLowerCase()]), n + r[0].length) : -1;
                            },
                            c: function (t, n, r) {
                                return _(t, e, n, r);
                            },
                            d: rn,
                            e: rn,
                            f: un,
                            H: on,
                            I: on,
                            j: an,
                            L: cn,
                            m: nn,
                            M: sn,
                            p: function (t, e, n) {
                                var r = c.exec(e.slice(n));
                                return r ? ((t.p = u[r[0].toLowerCase()]), n + r[0].length) : -1;
                            },
                            Q: dn,
                            s: gn,
                            S: ln,
                            u: Xe,
                            U: Ge,
                            V: Ke,
                            w: Ze,
                            W: Qe,
                            x: function (t, e, r) {
                                return _(t, n, e, r);
                            },
                            X: function (t, e, n) {
                                return _(t, r, e, n);
                            },
                            y: tn,
                            Y: Je,
                            Z: en,
                            "%": hn,
                        };
                    function $(t, e) {
                        return function (n) {
                            var r,
                                i,
                                a,
                                o = [],
                                s = -1,
                                l = 0,
                                c = t.length;
                            for (n instanceof Date || (n = new Date(+n)); ++s < c; )
                                37 === t.charCodeAt(s) &&
                                    (o.push(t.slice(l, s)),
                                    null != (i = ze[(r = t.charAt(++s))])
                                        ? (r = t.charAt(++s))
                                        : (i = "e" === r ? " " : "0"),
                                    (a = e[r]) && (r = a(n, i)),
                                    o.push(r),
                                    (l = s + 1));
                            return (o.push(t.slice(l, s)), o.join(""));
                        };
                    }
                    function A(t, e) {
                        return function (n) {
                            var r,
                                i,
                                a = Ie(1900);
                            if (_(a, t, (n += ""), 0) != n.length) return null;
                            if ("Q" in a) return new Date(a.Q);
                            if (("p" in a && (a.H = (a.H % 12) + 12 * a.p), "V" in a)) {
                                if (a.V < 1 || a.V > 53) return null;
                                ("w" in a || (a.w = 1),
                                    "Z" in a
                                        ? ((i = (r = De(Ie(a.y))).getUTCDay()),
                                          (r = i > 4 || 0 === i ? Ae.ceil(r) : Ae(r)),
                                          (r = xe.offset(r, 7 * (a.V - 1))),
                                          (a.y = r.getUTCFullYear()),
                                          (a.m = r.getUTCMonth()),
                                          (a.d = r.getUTCDate() + ((a.w + 6) % 7)))
                                        : ((i = (r = e(Ie(a.y))).getDay()),
                                          (r = i > 4 || 0 === i ? ae.ceil(r) : ae(r)),
                                          (r = ne.offset(r, 7 * (a.V - 1))),
                                          (a.y = r.getFullYear()),
                                          (a.m = r.getMonth()),
                                          (a.d = r.getDate() + ((a.w + 6) % 7))));
                            } else
                                ("W" in a || "U" in a) &&
                                    ("w" in a || (a.w = "u" in a ? a.u % 7 : "W" in a ? 1 : 0),
                                    (i = "Z" in a ? De(Ie(a.y)).getUTCDay() : e(Ie(a.y)).getDay()),
                                    (a.m = 0),
                                    (a.d =
                                        "W" in a
                                            ? ((a.w + 6) % 7) + 7 * a.W - ((i + 5) % 7)
                                            : a.w + 7 * a.U - ((i + 6) % 7)));
                            return "Z" in a ? ((a.H += (a.Z / 100) | 0), (a.M += a.Z % 100), De(a)) : e(a);
                        };
                    }
                    function _(t, e, n, r) {
                        for (var i, a, o = 0, s = e.length, l = n.length; o < s; ) {
                            if (r >= l) return -1;
                            if (37 === (i = e.charCodeAt(o++))) {
                                if (
                                    ((i = e.charAt(o++)), !(a = w[i in ze ? e.charAt(o++) : i]) || (r = a(t, n, r)) < 0)
                                )
                                    return -1;
                            } else if (i != n.charCodeAt(r++)) return -1;
                        }
                        return r;
                    }
                    return (
                        (b.x = $(n, b)),
                        (b.X = $(r, b)),
                        (b.c = $(e, b)),
                        (x.x = $(n, x)),
                        (x.X = $(r, x)),
                        (x.c = $(e, x)),
                        {
                            format: function (t) {
                                var e = $((t += ""), b);
                                return (
                                    (e.toString = function () {
                                        return t;
                                    }),
                                    e
                                );
                            },
                            parse: function (t) {
                                var e = A((t += ""), Pe);
                                return (
                                    (e.toString = function () {
                                        return t;
                                    }),
                                    e
                                );
                            },
                            utcFormat: function (t) {
                                var e = $((t += ""), x);
                                return (
                                    (e.toString = function () {
                                        return t;
                                    }),
                                    e
                                );
                            },
                            utcParse: function (t) {
                                var e = A(t, De);
                                return (
                                    (e.toString = function () {
                                        return t;
                                    }),
                                    e
                                );
                            },
                        }
                    );
                })({
                    dateTime: "%x, %X",
                    date: "%-m/%-d/%Y",
                    time: "%-I:%M:%S %p",
                    periods: ["AM", "PM"],
                    days: ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"],
                    shortDays: ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"],
                    months: [
                        "January",
                        "February",
                        "March",
                        "April",
                        "May",
                        "June",
                        "July",
                        "August",
                        "September",
                        "October",
                        "November",
                        "December",
                    ],
                    shortMonths: ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"],
                })),
                    (Be = Ee.format),
                    Ee.parse,
                    (Oe = Ee.utcFormat),
                    (Fe = Ee.utcParse));
                var Kn = "%Y-%m-%dT%H:%M:%S.%LZ";
                (Date.prototype.toISOString || Oe(Kn), +new Date("2000-01-01T00:00:00.000Z") || Fe(Kn));
                var Qn = 1e3,
                    Jn = 60 * Qn,
                    tr = 60 * Jn,
                    er = 24 * tr,
                    nr = 7 * er,
                    rr = 30 * er,
                    ir = 365 * er;
                function ar(t) {
                    return new Date(t);
                }
                function or(t) {
                    return t instanceof Date ? +t : +new Date(+t);
                }
                function sr(t, e, n, r, a, o, s, l, c) {
                    var u = J(W, W),
                        h = u.invert,
                        d = u.domain,
                        g = c(".%L"),
                        p = c(":%S"),
                        m = c("%I:%M"),
                        v = c("%I %p"),
                        y = c("%a %d"),
                        b = c("%b %d"),
                        x = c("%B"),
                        w = c("%Y"),
                        $ = [
                            [s, 1, Qn],
                            [s, 5, 5 * Qn],
                            [s, 15, 15 * Qn],
                            [s, 30, 30 * Qn],
                            [o, 1, Jn],
                            [o, 5, 5 * Jn],
                            [o, 15, 15 * Jn],
                            [o, 30, 30 * Jn],
                            [a, 1, tr],
                            [a, 3, 3 * tr],
                            [a, 6, 6 * tr],
                            [a, 12, 12 * tr],
                            [r, 1, er],
                            [r, 2, 2 * er],
                            [n, 1, nr],
                            [e, 1, rr],
                            [e, 3, 3 * rr],
                            [t, 1, ir],
                        ];
                    function A(i) {
                        return (
                            s(i) < i
                                ? g
                                : o(i) < i
                                  ? p
                                  : a(i) < i
                                    ? m
                                    : r(i) < i
                                      ? v
                                      : e(i) < i
                                        ? n(i) < i
                                            ? y
                                            : b
                                        : t(i) < i
                                          ? x
                                          : w
                        )(i);
                    }
                    function _(e, n, r, a) {
                        if ((null == e && (e = 10), "number" == typeof e)) {
                            var o = Math.abs(r - n) / e,
                                s = i(function (t) {
                                    return t[2];
                                }).right($, o);
                            s === $.length
                                ? ((a = f(n / ir, r / ir, e)), (e = t))
                                : s
                                  ? ((a = (s = $[o / $[s - 1][2] < $[s][2] / o ? s - 1 : s])[1]), (e = s[0]))
                                  : ((a = Math.max(f(n, r, e), 1)), (e = l));
                        }
                        return null == a ? e : e.every(a);
                    }
                    return (
                        (u.invert = function (t) {
                            return new Date(h(t));
                        }),
                        (u.domain = function (t) {
                            return arguments.length ? d(S.call(t, or)) : d().map(ar);
                        }),
                        (u.ticks = function (t, e) {
                            var n,
                                r = d(),
                                i = r[0],
                                a = r[r.length - 1],
                                o = a < i;
                            return (
                                o && ((n = i), (i = a), (a = n)),
                                (n = (n = _(t, i, a, e)) ? n.range(i, a + 1) : []),
                                o ? n.reverse() : n
                            );
                        }),
                        (u.tickFormat = function (t, e) {
                            return null == e ? A : c(e);
                        }),
                        (u.nice = function (t, e) {
                            var n = d();
                            return (t = _(t, n[0], n[n.length - 1], e)) ? d(vt(n, t)) : u;
                        }),
                        (u.copy = function () {
                            return K(u, sr(t, e, n, r, a, o, s, l, c));
                        }),
                        u
                    );
                }
                function lr() {
                    return v.apply(
                        sr(fe, de, ie, ne, te, Qt, Gt, jt, Be).domain([new Date(2e3, 0, 1), new Date(2e3, 0, 2)]),
                        arguments,
                    );
                }
                function cr() {
                    return v.apply(
                        sr(Le, Re, $e, xe, ye, me, Gt, jt, Oe).domain([Date.UTC(2e3, 0, 1), Date.UTC(2e3, 0, 2)]),
                        arguments,
                    );
                }
                function ur() {
                    var t,
                        e,
                        n,
                        r,
                        i,
                        a = 0,
                        o = 1,
                        s = W,
                        l = !1;
                    function c(e) {
                        return isNaN((e = +e))
                            ? i
                            : s(0 === n ? 0.5 : ((e = (r(e) - t) * n), l ? Math.max(0, Math.min(1, e)) : e));
                    }
                    return (
                        (c.domain = function (i) {
                            return arguments.length
                                ? ((t = r((a = +i[0]))), (e = r((o = +i[1]))), (n = t === e ? 0 : 1 / (e - t)), c)
                                : [a, o];
                        }),
                        (c.clamp = function (t) {
                            return arguments.length ? ((l = !!t), c) : l;
                        }),
                        (c.interpolator = function (t) {
                            return arguments.length ? ((s = t), c) : s;
                        }),
                        (c.unknown = function (t) {
                            return arguments.length ? ((i = t), c) : i;
                        }),
                        function (i) {
                            return ((r = i), (t = i(a)), (e = i(o)), (n = t === e ? 0 : 1 / (e - t)), c);
                        }
                    );
                }
                function hr(t, e) {
                    return e.domain(t.domain()).interpolator(t.interpolator()).clamp(t.clamp()).unknown(t.unknown());
                }
                function dr() {
                    var t = ft(ur()(W));
                    return (
                        (t.copy = function () {
                            return hr(t, dr());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function gr() {
                    var t = _t(ur()).domain([1, 10]);
                    return (
                        (t.copy = function () {
                            return hr(t, gr()).base(t.base());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function fr() {
                    var t = Tt(ur());
                    return (
                        (t.copy = function () {
                            return hr(t, fr()).constant(t.constant());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function pr() {
                    var t = Pt(ur());
                    return (
                        (t.copy = function () {
                            return hr(t, pr()).exponent(t.exponent());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function mr() {
                    return pr.apply(null, arguments).exponent(0.5);
                }
                function vr() {
                    var t = [],
                        e = W;
                    function n(n) {
                        if (!isNaN((n = +n))) return e((s(t, n) - 1) / (t.length - 1));
                    }
                    return (
                        (n.domain = function (e) {
                            if (!arguments.length) return t.slice();
                            t = [];
                            for (var i, a = 0, o = e.length; a < o; ++a)
                                null == (i = e[a]) || isNaN((i = +i)) || t.push(i);
                            return (t.sort(r), n);
                        }),
                        (n.interpolator = function (t) {
                            return arguments.length ? ((e = t), n) : e;
                        }),
                        (n.copy = function () {
                            return vr(e).domain(t);
                        }),
                        y.apply(n, arguments)
                    );
                }
                function yr() {
                    var t,
                        e,
                        n,
                        r,
                        i,
                        a,
                        o,
                        s = 0,
                        l = 0.5,
                        c = 1,
                        u = W,
                        h = !1;
                    function d(t) {
                        return isNaN((t = +t))
                            ? o
                            : ((t = 0.5 + ((t = +a(t)) - e) * (t < e ? r : i)), u(h ? Math.max(0, Math.min(1, t)) : t));
                    }
                    return (
                        (d.domain = function (o) {
                            return arguments.length
                                ? ((t = a((s = +o[0]))),
                                  (e = a((l = +o[1]))),
                                  (n = a((c = +o[2]))),
                                  (r = t === e ? 0 : 0.5 / (e - t)),
                                  (i = e === n ? 0 : 0.5 / (n - e)),
                                  d)
                                : [s, l, c];
                        }),
                        (d.clamp = function (t) {
                            return arguments.length ? ((h = !!t), d) : h;
                        }),
                        (d.interpolator = function (t) {
                            return arguments.length ? ((u = t), d) : u;
                        }),
                        (d.unknown = function (t) {
                            return arguments.length ? ((o = t), d) : o;
                        }),
                        function (o) {
                            return (
                                (a = o),
                                (t = o(s)),
                                (e = o(l)),
                                (n = o(c)),
                                (r = t === e ? 0 : 0.5 / (e - t)),
                                (i = e === n ? 0 : 0.5 / (n - e)),
                                d
                            );
                        }
                    );
                }
                function br() {
                    var t = ft(yr()(W));
                    return (
                        (t.copy = function () {
                            return hr(t, br());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function xr() {
                    var t = _t(yr()).domain([0.1, 1, 10]);
                    return (
                        (t.copy = function () {
                            return hr(t, xr()).base(t.base());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function wr() {
                    var t = Tt(yr());
                    return (
                        (t.copy = function () {
                            return hr(t, wr()).constant(t.constant());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function $r() {
                    var t = Pt(yr());
                    return (
                        (t.copy = function () {
                            return hr(t, $r()).exponent(t.exponent());
                        }),
                        y.apply(t, arguments)
                    );
                }
                function Ar() {
                    return $r.apply(null, arguments).exponent(0.5);
                }
            },
            677: (t, e, n) => {
                "use strict";
                (n.r(e),
                    n.d(e, {
                        clientPoint: () => $t,
                        create: () => vt,
                        creator: () => l,
                        customEvent: () => lt,
                        event: () => rt,
                        local: () => bt,
                        matcher: () => m,
                        mouse: () => At,
                        namespace: () => a,
                        namespaces: () => i,
                        select: () => mt,
                        selectAll: () => _t,
                        selection: () => pt,
                        selector: () => u,
                        selectorAll: () => d,
                        style: () => L,
                        touch: () => Ct,
                        touches: () => St,
                        window: () => T,
                    }));
                var r = "http://www.w3.org/1999/xhtml";
                const i = {
                    svg: "http://www.w3.org/2000/svg",
                    xhtml: r,
                    xlink: "http://www.w3.org/1999/xlink",
                    xml: "http://www.w3.org/XML/1998/namespace",
                    xmlns: "http://www.w3.org/2000/xmlns/",
                };
                function a(t) {
                    var e = (t += ""),
                        n = e.indexOf(":");
                    return (
                        n >= 0 && "xmlns" !== (e = t.slice(0, n)) && (t = t.slice(n + 1)),
                        i.hasOwnProperty(e) ? { space: i[e], local: t } : t
                    );
                }
                function o(t) {
                    return function () {
                        var e = this.ownerDocument,
                            n = this.namespaceURI;
                        return n === r && e.documentElement.namespaceURI === r
                            ? e.createElement(t)
                            : e.createElementNS(n, t);
                    };
                }
                function s(t) {
                    return function () {
                        return this.ownerDocument.createElementNS(t.space, t.local);
                    };
                }
                function l(t) {
                    var e = a(t);
                    return (e.local ? s : o)(e);
                }
                function c() {}
                function u(t) {
                    return null == t
                        ? c
                        : function () {
                              return this.querySelector(t);
                          };
                }
                function h() {
                    return [];
                }
                function d(t) {
                    return null == t
                        ? h
                        : function () {
                              return this.querySelectorAll(t);
                          };
                }
                var g = function (t) {
                    return function () {
                        return this.matches(t);
                    };
                };
                if ("undefined" != typeof document) {
                    var f = document.documentElement;
                    if (!f.matches) {
                        var p =
                            f.webkitMatchesSelector ||
                            f.msMatchesSelector ||
                            f.mozMatchesSelector ||
                            f.oMatchesSelector;
                        g = function (t) {
                            return function () {
                                return p.call(this, t);
                            };
                        };
                    }
                }
                const m = g;
                function v(t) {
                    return new Array(t.length);
                }
                function y(t, e) {
                    ((this.ownerDocument = t.ownerDocument),
                        (this.namespaceURI = t.namespaceURI),
                        (this._next = null),
                        (this._parent = t),
                        (this.__data__ = e));
                }
                y.prototype = {
                    constructor: y,
                    appendChild: function (t) {
                        return this._parent.insertBefore(t, this._next);
                    },
                    insertBefore: function (t, e) {
                        return this._parent.insertBefore(t, e);
                    },
                    querySelector: function (t) {
                        return this._parent.querySelector(t);
                    },
                    querySelectorAll: function (t) {
                        return this._parent.querySelectorAll(t);
                    },
                };
                function b(t, e, n, r, i, a) {
                    for (var o, s = 0, l = e.length, c = a.length; s < c; ++s)
                        (o = e[s]) ? ((o.__data__ = a[s]), (r[s] = o)) : (n[s] = new y(t, a[s]));
                    for (; s < l; ++s) (o = e[s]) && (i[s] = o);
                }
                function x(t, e, n, r, i, a, o) {
                    var s,
                        l,
                        c,
                        u = {},
                        h = e.length,
                        d = a.length,
                        g = new Array(h);
                    for (s = 0; s < h; ++s)
                        (l = e[s]) &&
                            ((g[s] = c = "$" + o.call(l, l.__data__, s, e)), c in u ? (i[s] = l) : (u[c] = l));
                    for (s = 0; s < d; ++s)
                        (l = u[(c = "$" + o.call(t, a[s], s, a))])
                            ? ((r[s] = l), (l.__data__ = a[s]), (u[c] = null))
                            : (n[s] = new y(t, a[s]));
                    for (s = 0; s < h; ++s) (l = e[s]) && u[g[s]] === l && (i[s] = l);
                }
                function w(t, e) {
                    return t < e ? -1 : t > e ? 1 : t >= e ? 0 : NaN;
                }
                function $(t) {
                    return function () {
                        this.removeAttribute(t);
                    };
                }
                function A(t) {
                    return function () {
                        this.removeAttributeNS(t.space, t.local);
                    };
                }
                function _(t, e) {
                    return function () {
                        this.setAttribute(t, e);
                    };
                }
                function C(t, e) {
                    return function () {
                        this.setAttributeNS(t.space, t.local, e);
                    };
                }
                function S(t, e) {
                    return function () {
                        var n = e.apply(this, arguments);
                        null == n ? this.removeAttribute(t) : this.setAttribute(t, n);
                    };
                }
                function M(t, e) {
                    return function () {
                        var n = e.apply(this, arguments);
                        null == n ? this.removeAttributeNS(t.space, t.local) : this.setAttributeNS(t.space, t.local, n);
                    };
                }
                function T(t) {
                    return (t.ownerDocument && t.ownerDocument.defaultView) || (t.document && t) || t.defaultView;
                }
                function k(t) {
                    return function () {
                        this.style.removeProperty(t);
                    };
                }
                function R(t, e, n) {
                    return function () {
                        this.style.setProperty(t, e, n);
                    };
                }
                function N(t, e, n) {
                    return function () {
                        var r = e.apply(this, arguments);
                        null == r ? this.style.removeProperty(t) : this.style.setProperty(t, r, n);
                    };
                }
                function L(t, e) {
                    return t.style.getPropertyValue(e) || T(t).getComputedStyle(t, null).getPropertyValue(e);
                }
                function P(t) {
                    return function () {
                        delete this[t];
                    };
                }
                function D(t, e) {
                    return function () {
                        this[t] = e;
                    };
                }
                function I(t, e) {
                    return function () {
                        var n = e.apply(this, arguments);
                        null == n ? delete this[t] : (this[t] = n);
                    };
                }
                function E(t) {
                    return t.trim().split(/^|\s+/);
                }
                function B(t) {
                    return t.classList || new O(t);
                }
                function O(t) {
                    ((this._node = t), (this._names = E(t.getAttribute("class") || "")));
                }
                function F(t, e) {
                    for (var n = B(t), r = -1, i = e.length; ++r < i; ) n.add(e[r]);
                }
                function z(t, e) {
                    for (var n = B(t), r = -1, i = e.length; ++r < i; ) n.remove(e[r]);
                }
                function q(t) {
                    return function () {
                        F(this, t);
                    };
                }
                function H(t) {
                    return function () {
                        z(this, t);
                    };
                }
                function j(t, e) {
                    return function () {
                        (e.apply(this, arguments) ? F : z)(this, t);
                    };
                }
                function U() {
                    this.textContent = "";
                }
                function V(t) {
                    return function () {
                        this.textContent = t;
                    };
                }
                function W(t) {
                    return function () {
                        var e = t.apply(this, arguments);
                        this.textContent = null == e ? "" : e;
                    };
                }
                function Y() {
                    this.innerHTML = "";
                }
                function Z(t) {
                    return function () {
                        this.innerHTML = t;
                    };
                }
                function X(t) {
                    return function () {
                        var e = t.apply(this, arguments);
                        this.innerHTML = null == e ? "" : e;
                    };
                }
                function G() {
                    this.nextSibling && this.parentNode.appendChild(this);
                }
                function K() {
                    this.previousSibling && this.parentNode.insertBefore(this, this.parentNode.firstChild);
                }
                function Q() {
                    return null;
                }
                function J() {
                    var t = this.parentNode;
                    t && t.removeChild(this);
                }
                function tt() {
                    return this.parentNode.insertBefore(this.cloneNode(!1), this.nextSibling);
                }
                function et() {
                    return this.parentNode.insertBefore(this.cloneNode(!0), this.nextSibling);
                }
                O.prototype = {
                    add: function (t) {
                        this._names.indexOf(t) < 0 &&
                            (this._names.push(t), this._node.setAttribute("class", this._names.join(" ")));
                    },
                    remove: function (t) {
                        var e = this._names.indexOf(t);
                        e >= 0 && (this._names.splice(e, 1), this._node.setAttribute("class", this._names.join(" ")));
                    },
                    contains: function (t) {
                        return this._names.indexOf(t) >= 0;
                    },
                };
                var nt = {},
                    rt = null;
                function it(t, e, n) {
                    return (
                        (t = at(t, e, n)),
                        function (e) {
                            var n = e.relatedTarget;
                            (n && (n === this || 8 & n.compareDocumentPosition(this))) || t.call(this, e);
                        }
                    );
                }
                function at(t, e, n) {
                    return function (r) {
                        var i = rt;
                        rt = r;
                        try {
                            t.call(this, this.__data__, e, n);
                        } finally {
                            rt = i;
                        }
                    };
                }
                function ot(t) {
                    return function () {
                        var e = this.__on;
                        if (e) {
                            for (var n, r = 0, i = -1, a = e.length; r < a; ++r)
                                ((n = e[r]),
                                    (t.type && n.type !== t.type) || n.name !== t.name
                                        ? (e[++i] = n)
                                        : this.removeEventListener(n.type, n.listener, n.capture));
                            ++i ? (e.length = i) : delete this.__on;
                        }
                    };
                }
                function st(t, e, n) {
                    var r = nt.hasOwnProperty(t.type) ? it : at;
                    return function (i, a, o) {
                        var s,
                            l = this.__on,
                            c = r(e, a, o);
                        if (l)
                            for (var u = 0, h = l.length; u < h; ++u)
                                if ((s = l[u]).type === t.type && s.name === t.name)
                                    return (
                                        this.removeEventListener(s.type, s.listener, s.capture),
                                        this.addEventListener(s.type, (s.listener = c), (s.capture = n)),
                                        void (s.value = e)
                                    );
                        (this.addEventListener(t.type, c, n),
                            (s = { type: t.type, name: t.name, value: e, listener: c, capture: n }),
                            l ? l.push(s) : (this.__on = [s]));
                    };
                }
                function lt(t, e, n, r) {
                    var i = rt;
                    ((t.sourceEvent = rt), (rt = t));
                    try {
                        return e.apply(n, r);
                    } finally {
                        rt = i;
                    }
                }
                function ct(t, e, n) {
                    var r = T(t),
                        i = r.CustomEvent;
                    ("function" == typeof i
                        ? (i = new i(e, n))
                        : ((i = r.document.createEvent("Event")),
                          n
                              ? (i.initEvent(e, n.bubbles, n.cancelable), (i.detail = n.detail))
                              : i.initEvent(e, !1, !1)),
                        t.dispatchEvent(i));
                }
                function ut(t, e) {
                    return function () {
                        return ct(this, t, e);
                    };
                }
                function ht(t, e) {
                    return function () {
                        return ct(this, t, e.apply(this, arguments));
                    };
                }
                "undefined" != typeof document &&
                    ("onmouseenter" in document.documentElement ||
                        (nt = { mouseenter: "mouseover", mouseleave: "mouseout" }));
                var dt = [null];
                function gt(t, e) {
                    ((this._groups = t), (this._parents = e));
                }
                function ft() {
                    return new gt([[document.documentElement]], dt);
                }
                gt.prototype = ft.prototype = {
                    constructor: gt,
                    select: function (t) {
                        "function" != typeof t && (t = u(t));
                        for (var e = this._groups, n = e.length, r = new Array(n), i = 0; i < n; ++i)
                            for (var a, o, s = e[i], l = s.length, c = (r[i] = new Array(l)), h = 0; h < l; ++h)
                                (a = s[h]) &&
                                    (o = t.call(a, a.__data__, h, s)) &&
                                    ("__data__" in a && (o.__data__ = a.__data__), (c[h] = o));
                        return new gt(r, this._parents);
                    },
                    selectAll: function (t) {
                        "function" != typeof t && (t = d(t));
                        for (var e = this._groups, n = e.length, r = [], i = [], a = 0; a < n; ++a)
                            for (var o, s = e[a], l = s.length, c = 0; c < l; ++c)
                                (o = s[c]) && (r.push(t.call(o, o.__data__, c, s)), i.push(o));
                        return new gt(r, i);
                    },
                    filter: function (t) {
                        "function" != typeof t && (t = m(t));
                        for (var e = this._groups, n = e.length, r = new Array(n), i = 0; i < n; ++i)
                            for (var a, o = e[i], s = o.length, l = (r[i] = []), c = 0; c < s; ++c)
                                (a = o[c]) && t.call(a, a.__data__, c, o) && l.push(a);
                        return new gt(r, this._parents);
                    },
                    data: function (t, e) {
                        if (!t)
                            return (
                                (f = new Array(this.size())),
                                (u = -1),
                                this.each(function (t) {
                                    f[++u] = t;
                                }),
                                f
                            );
                        var n,
                            r = e ? x : b,
                            i = this._parents,
                            a = this._groups;
                        "function" != typeof t &&
                            ((n = t),
                            (t = function () {
                                return n;
                            }));
                        for (
                            var o = a.length, s = new Array(o), l = new Array(o), c = new Array(o), u = 0;
                            u < o;
                            ++u
                        ) {
                            var h = i[u],
                                d = a[u],
                                g = d.length,
                                f = t.call(h, h && h.__data__, u, i),
                                p = f.length,
                                m = (l[u] = new Array(p)),
                                v = (s[u] = new Array(p));
                            r(h, d, m, v, (c[u] = new Array(g)), f, e);
                            for (var y, w, $ = 0, A = 0; $ < p; ++$)
                                if ((y = m[$])) {
                                    for ($ >= A && (A = $ + 1); !(w = v[A]) && ++A < p; );
                                    y._next = w || null;
                                }
                        }
                        return (((s = new gt(s, i))._enter = l), (s._exit = c), s);
                    },
                    enter: function () {
                        return new gt(this._enter || this._groups.map(v), this._parents);
                    },
                    exit: function () {
                        return new gt(this._exit || this._groups.map(v), this._parents);
                    },
                    merge: function (t) {
                        for (
                            var e = this._groups,
                                n = t._groups,
                                r = e.length,
                                i = n.length,
                                a = Math.min(r, i),
                                o = new Array(r),
                                s = 0;
                            s < a;
                            ++s
                        )
                            for (var l, c = e[s], u = n[s], h = c.length, d = (o[s] = new Array(h)), g = 0; g < h; ++g)
                                (l = c[g] || u[g]) && (d[g] = l);
                        for (; s < r; ++s) o[s] = e[s];
                        return new gt(o, this._parents);
                    },
                    order: function () {
                        for (var t = this._groups, e = -1, n = t.length; ++e < n; )
                            for (var r, i = t[e], a = i.length - 1, o = i[a]; --a >= 0; )
                                (r = i[a]) && (o && o !== r.nextSibling && o.parentNode.insertBefore(r, o), (o = r));
                        return this;
                    },
                    sort: function (t) {
                        function e(e, n) {
                            return e && n ? t(e.__data__, n.__data__) : !e - !n;
                        }
                        t || (t = w);
                        for (var n = this._groups, r = n.length, i = new Array(r), a = 0; a < r; ++a) {
                            for (var o, s = n[a], l = s.length, c = (i[a] = new Array(l)), u = 0; u < l; ++u)
                                (o = s[u]) && (c[u] = o);
                            c.sort(e);
                        }
                        return new gt(i, this._parents).order();
                    },
                    call: function () {
                        var t = arguments[0];
                        return ((arguments[0] = this), t.apply(null, arguments), this);
                    },
                    nodes: function () {
                        var t = new Array(this.size()),
                            e = -1;
                        return (
                            this.each(function () {
                                t[++e] = this;
                            }),
                            t
                        );
                    },
                    node: function () {
                        for (var t = this._groups, e = 0, n = t.length; e < n; ++e)
                            for (var r = t[e], i = 0, a = r.length; i < a; ++i) {
                                var o = r[i];
                                if (o) return o;
                            }
                        return null;
                    },
                    size: function () {
                        var t = 0;
                        return (
                            this.each(function () {
                                ++t;
                            }),
                            t
                        );
                    },
                    empty: function () {
                        return !this.node();
                    },
                    each: function (t) {
                        for (var e = this._groups, n = 0, r = e.length; n < r; ++n)
                            for (var i, a = e[n], o = 0, s = a.length; o < s; ++o)
                                (i = a[o]) && t.call(i, i.__data__, o, a);
                        return this;
                    },
                    attr: function (t, e) {
                        var n = a(t);
                        if (arguments.length < 2) {
                            var r = this.node();
                            return n.local ? r.getAttributeNS(n.space, n.local) : r.getAttribute(n);
                        }
                        return this.each(
                            (null == e
                                ? n.local
                                    ? A
                                    : $
                                : "function" == typeof e
                                  ? n.local
                                      ? M
                                      : S
                                  : n.local
                                    ? C
                                    : _)(n, e),
                        );
                    },
                    style: function (t, e, n) {
                        return arguments.length > 1
                            ? this.each((null == e ? k : "function" == typeof e ? N : R)(t, e, null == n ? "" : n))
                            : L(this.node(), t);
                    },
                    property: function (t, e) {
                        return arguments.length > 1
                            ? this.each((null == e ? P : "function" == typeof e ? I : D)(t, e))
                            : this.node()[t];
                    },
                    classed: function (t, e) {
                        var n = E(t + "");
                        if (arguments.length < 2) {
                            for (var r = B(this.node()), i = -1, a = n.length; ++i < a; )
                                if (!r.contains(n[i])) return !1;
                            return !0;
                        }
                        return this.each(("function" == typeof e ? j : e ? q : H)(n, e));
                    },
                    text: function (t) {
                        return arguments.length
                            ? this.each(null == t ? U : ("function" == typeof t ? W : V)(t))
                            : this.node().textContent;
                    },
                    html: function (t) {
                        return arguments.length
                            ? this.each(null == t ? Y : ("function" == typeof t ? X : Z)(t))
                            : this.node().innerHTML;
                    },
                    raise: function () {
                        return this.each(G);
                    },
                    lower: function () {
                        return this.each(K);
                    },
                    append: function (t) {
                        var e = "function" == typeof t ? t : l(t);
                        return this.select(function () {
                            return this.appendChild(e.apply(this, arguments));
                        });
                    },
                    insert: function (t, e) {
                        var n = "function" == typeof t ? t : l(t),
                            r = null == e ? Q : "function" == typeof e ? e : u(e);
                        return this.select(function () {
                            return this.insertBefore(n.apply(this, arguments), r.apply(this, arguments) || null);
                        });
                    },
                    remove: function () {
                        return this.each(J);
                    },
                    clone: function (t) {
                        return this.select(t ? et : tt);
                    },
                    datum: function (t) {
                        return arguments.length ? this.property("__data__", t) : this.node().__data__;
                    },
                    on: function (t, e, n) {
                        var r,
                            i,
                            a = (function (t) {
                                return t
                                    .trim()
                                    .split(/^|\s+/)
                                    .map(function (t) {
                                        var e = "",
                                            n = t.indexOf(".");
                                        return (
                                            n >= 0 && ((e = t.slice(n + 1)), (t = t.slice(0, n))),
                                            { type: t, name: e }
                                        );
                                    });
                            })(t + ""),
                            o = a.length;
                        if (!(arguments.length < 2)) {
                            for (s = e ? st : ot, null == n && (n = !1), r = 0; r < o; ++r) this.each(s(a[r], e, n));
                            return this;
                        }
                        var s = this.node().__on;
                        if (s)
                            for (var l, c = 0, u = s.length; c < u; ++c)
                                for (r = 0, l = s[c]; r < o; ++r)
                                    if ((i = a[r]).type === l.type && i.name === l.name) return l.value;
                    },
                    dispatch: function (t, e) {
                        return this.each(("function" == typeof e ? ht : ut)(t, e));
                    },
                };
                const pt = ft;
                function mt(t) {
                    return "string" == typeof t
                        ? new gt([[document.querySelector(t)]], [document.documentElement])
                        : new gt([[t]], dt);
                }
                function vt(t) {
                    return mt(l(t).call(document.documentElement));
                }
                var yt = 0;
                function bt() {
                    return new xt();
                }
                function xt() {
                    this._ = "@" + (++yt).toString(36);
                }
                function wt() {
                    for (var t, e = rt; (t = e.sourceEvent); ) e = t;
                    return e;
                }
                function $t(t, e) {
                    var n = t.ownerSVGElement || t;
                    if (n.createSVGPoint) {
                        var r = n.createSVGPoint();
                        return (
                            (r.x = e.clientX),
                            (r.y = e.clientY),
                            [(r = r.matrixTransform(t.getScreenCTM().inverse())).x, r.y]
                        );
                    }
                    var i = t.getBoundingClientRect();
                    return [e.clientX - i.left - t.clientLeft, e.clientY - i.top - t.clientTop];
                }
                function At(t) {
                    var e = wt();
                    return (e.changedTouches && (e = e.changedTouches[0]), $t(t, e));
                }
                function _t(t) {
                    return "string" == typeof t
                        ? new gt([document.querySelectorAll(t)], [document.documentElement])
                        : new gt([null == t ? [] : t], dt);
                }
                function Ct(t, e, n) {
                    arguments.length < 3 && ((n = e), (e = wt().changedTouches));
                    for (var r, i = 0, a = e ? e.length : 0; i < a; ++i)
                        if ((r = e[i]).identifier === n) return $t(t, r);
                    return null;
                }
                function St(t, e) {
                    null == e && (e = wt().touches);
                    for (var n = 0, r = e ? e.length : 0, i = new Array(r); n < r; ++n) i[n] = $t(t, e[n]);
                    return i;
                }
                xt.prototype = bt.prototype = {
                    constructor: xt,
                    get: function (t) {
                        for (var e = this._; !(e in t); ) if (!(t = t.parentNode)) return;
                        return t[e];
                    },
                    set: function (t, e) {
                        return (t[this._] = e);
                    },
                    remove: function (t) {
                        return this._ in t && delete t[this._];
                    },
                    toString: function () {
                        return this._;
                    },
                };
            },
            473: (t, e, n) => {
                "use strict";
                (n.r(e), n.d(e, { active: () => Ct, interrupt: () => E, transition: () => vt }));
                var r,
                    i,
                    a = n(677),
                    o = n(594),
                    s = 0,
                    l = 0,
                    c = 0,
                    u = 1e3,
                    h = 0,
                    d = 0,
                    g = 0,
                    f = "object" == typeof performance && performance.now ? performance : Date,
                    p =
                        "object" == typeof window && window.requestAnimationFrame
                            ? window.requestAnimationFrame.bind(window)
                            : function (t) {
                                  setTimeout(t, 17);
                              };
                function m() {
                    return d || (p(v), (d = f.now() + g));
                }
                function v() {
                    d = 0;
                }
                function y() {
                    this._call = this._time = this._next = null;
                }
                function b(t, e, n) {
                    var r = new y();
                    return (r.restart(t, e, n), r);
                }
                function x() {
                    ((d = (h = f.now()) + g), (s = l = 0));
                    try {
                        !(function () {
                            (m(), ++s);
                            for (var t, e = r; e; ) ((t = d - e._time) >= 0 && e._call.call(null, t), (e = e._next));
                            --s;
                        })();
                    } finally {
                        ((s = 0),
                            (function () {
                                for (var t, e, n = r, a = 1 / 0; n; )
                                    n._call
                                        ? (a > n._time && (a = n._time), (t = n), (n = n._next))
                                        : ((e = n._next), (n._next = null), (n = t ? (t._next = e) : (r = e)));
                                ((i = t), $(a));
                            })(),
                            (d = 0));
                    }
                }
                function w() {
                    var t = f.now(),
                        e = t - h;
                    e > u && ((g -= e), (h = t));
                }
                function $(t) {
                    s ||
                        (l && (l = clearTimeout(l)),
                        t - d > 24
                            ? (t < 1 / 0 && (l = setTimeout(x, t - f.now() - g)), c && (c = clearInterval(c)))
                            : (c || ((h = f.now()), (c = setInterval(w, u))), (s = 1), p(x)));
                }
                function A(t, e, n) {
                    var r = new y();
                    return (
                        (e = null == e ? 0 : +e),
                        r.restart(
                            function (n) {
                                (r.stop(), t(n + e));
                            },
                            e,
                            n,
                        ),
                        r
                    );
                }
                y.prototype = b.prototype = {
                    constructor: y,
                    restart: function (t, e, n) {
                        if ("function" != typeof t) throw new TypeError("callback is not a function");
                        ((n = (null == n ? m() : +n) + (null == e ? 0 : +e)),
                            this._next || i === this || (i ? (i._next = this) : (r = this), (i = this)),
                            (this._call = t),
                            (this._time = n),
                            $());
                    },
                    stop: function () {
                        this._call && ((this._call = null), (this._time = 1 / 0), $());
                    },
                };
                var _ = (0, o.W)("start", "end", "cancel", "interrupt"),
                    C = [],
                    S = 0,
                    M = 1,
                    T = 2,
                    k = 3,
                    R = 5,
                    N = 6;
                function L(t, e, n, r, i, a) {
                    var o = t.__transition;
                    if (o) {
                        if (n in o) return;
                    } else t.__transition = {};
                    !(function (t, e, n) {
                        var r,
                            i = t.__transition;
                        function a(l) {
                            var c, u, h, d;
                            if (n.state !== M) return s();
                            for (c in i)
                                if ((d = i[c]).name === n.name) {
                                    if (d.state === k) return A(a);
                                    4 === d.state
                                        ? ((d.state = N),
                                          d.timer.stop(),
                                          d.on.call("interrupt", t, t.__data__, d.index, d.group),
                                          delete i[c])
                                        : +c < e &&
                                          ((d.state = N),
                                          d.timer.stop(),
                                          d.on.call("cancel", t, t.__data__, d.index, d.group),
                                          delete i[c]);
                                }
                            if (
                                (A(function () {
                                    n.state === k && ((n.state = 4), n.timer.restart(o, n.delay, n.time), o(l));
                                }),
                                (n.state = T),
                                n.on.call("start", t, t.__data__, n.index, n.group),
                                n.state === T)
                            ) {
                                for (n.state = k, r = new Array((h = n.tween.length)), c = 0, u = -1; c < h; ++c)
                                    (d = n.tween[c].value.call(t, t.__data__, n.index, n.group)) && (r[++u] = d);
                                r.length = u + 1;
                            }
                        }
                        function o(e) {
                            for (
                                var i =
                                        e < n.duration
                                            ? n.ease.call(null, e / n.duration)
                                            : (n.timer.restart(s), (n.state = R), 1),
                                    a = -1,
                                    o = r.length;
                                ++a < o;
                            )
                                r[a].call(t, i);
                            n.state === R && (n.on.call("end", t, t.__data__, n.index, n.group), s());
                        }
                        function s() {
                            for (var r in ((n.state = N), n.timer.stop(), delete i[e], i)) return;
                            delete t.__transition;
                        }
                        ((i[e] = n),
                            (n.timer = b(
                                function (t) {
                                    ((n.state = M),
                                        n.timer.restart(a, n.delay, n.time),
                                        n.delay <= t && a(t - n.delay));
                                },
                                0,
                                n.time,
                            )));
                    })(t, n, {
                        name: e,
                        index: r,
                        group: i,
                        on: _,
                        tween: C,
                        time: a.time,
                        delay: a.delay,
                        duration: a.duration,
                        ease: a.ease,
                        timer: null,
                        state: S,
                    });
                }
                function P(t, e) {
                    var n = I(t, e);
                    if (n.state > S) throw new Error("too late; already scheduled");
                    return n;
                }
                function D(t, e) {
                    var n = I(t, e);
                    if (n.state > k) throw new Error("too late; already running");
                    return n;
                }
                function I(t, e) {
                    var n = t.__transition;
                    if (!n || !(n = n[e])) throw new Error("transition not found");
                    return n;
                }
                function E(t, e) {
                    var n,
                        r,
                        i,
                        a = t.__transition,
                        o = !0;
                    if (a) {
                        for (i in ((e = null == e ? null : e + ""), a))
                            (n = a[i]).name === e
                                ? ((r = n.state > T && n.state < R),
                                  (n.state = N),
                                  n.timer.stop(),
                                  n.on.call(r ? "interrupt" : "cancel", t, t.__data__, n.index, n.group),
                                  delete a[i])
                                : (o = !1);
                        o && delete t.__transition;
                    }
                }
                var B,
                    O,
                    F,
                    z,
                    q = n(626),
                    H = 180 / Math.PI,
                    j = { translateX: 0, translateY: 0, rotate: 0, skewX: 0, scaleX: 1, scaleY: 1 };
                function U(t, e, n, r, i, a) {
                    var o, s, l;
                    return (
                        (o = Math.sqrt(t * t + e * e)) && ((t /= o), (e /= o)),
                        (l = t * n + e * r) && ((n -= t * l), (r -= e * l)),
                        (s = Math.sqrt(n * n + r * r)) && ((n /= s), (r /= s), (l /= s)),
                        t * r < e * n && ((t = -t), (e = -e), (l = -l), (o = -o)),
                        {
                            translateX: i,
                            translateY: a,
                            rotate: Math.atan2(e, t) * H,
                            skewX: Math.atan(l) * H,
                            scaleX: o,
                            scaleY: s,
                        }
                    );
                }
                function V(t, e, n, r) {
                    function i(t) {
                        return t.length ? t.pop() + " " : "";
                    }
                    return function (a, o) {
                        var s = [],
                            l = [];
                        return (
                            (a = t(a)),
                            (o = t(o)),
                            (function (t, r, i, a, o, s) {
                                if (t !== i || r !== a) {
                                    var l = o.push("translate(", null, e, null, n);
                                    s.push({ i: l - 4, x: (0, q.Z)(t, i) }, { i: l - 2, x: (0, q.Z)(r, a) });
                                } else (i || a) && o.push("translate(" + i + e + a + n);
                            })(a.translateX, a.translateY, o.translateX, o.translateY, s, l),
                            (function (t, e, n, a) {
                                t !== e
                                    ? (t - e > 180 ? (e += 360) : e - t > 180 && (t += 360),
                                      a.push({ i: n.push(i(n) + "rotate(", null, r) - 2, x: (0, q.Z)(t, e) }))
                                    : e && n.push(i(n) + "rotate(" + e + r);
                            })(a.rotate, o.rotate, s, l),
                            (function (t, e, n, a) {
                                t !== e
                                    ? a.push({ i: n.push(i(n) + "skewX(", null, r) - 2, x: (0, q.Z)(t, e) })
                                    : e && n.push(i(n) + "skewX(" + e + r);
                            })(a.skewX, o.skewX, s, l),
                            (function (t, e, n, r, a, o) {
                                if (t !== n || e !== r) {
                                    var s = a.push(i(a) + "scale(", null, ",", null, ")");
                                    o.push({ i: s - 4, x: (0, q.Z)(t, n) }, { i: s - 2, x: (0, q.Z)(e, r) });
                                } else (1 === n && 1 === r) || a.push(i(a) + "scale(" + n + "," + r + ")");
                            })(a.scaleX, a.scaleY, o.scaleX, o.scaleY, s, l),
                            (a = o = null),
                            function (t) {
                                for (var e, n = -1, r = l.length; ++n < r; ) s[(e = l[n]).i] = e.x(t);
                                return s.join("");
                            }
                        );
                    };
                }
                var W = V(
                        function (t) {
                            return "none" === t
                                ? j
                                : (B ||
                                      ((B = document.createElement("DIV")),
                                      (O = document.documentElement),
                                      (F = document.defaultView)),
                                  (B.style.transform = t),
                                  (t = F.getComputedStyle(O.appendChild(B), null).getPropertyValue("transform")),
                                  O.removeChild(B),
                                  U(+(t = t.slice(7, -1).split(","))[0], +t[1], +t[2], +t[3], +t[4], +t[5]));
                        },
                        "px, ",
                        "px)",
                        "deg)",
                    ),
                    Y = V(
                        function (t) {
                            return null == t
                                ? j
                                : (z || (z = document.createElementNS("http://www.w3.org/2000/svg", "g")),
                                  z.setAttribute("transform", t),
                                  (t = z.transform.baseVal.consolidate())
                                      ? U((t = t.matrix).a, t.b, t.c, t.d, t.e, t.f)
                                      : j);
                        },
                        ", ",
                        ")",
                        ")",
                    );
                function Z(t, e) {
                    var n, r;
                    return function () {
                        var i = D(this, t),
                            a = i.tween;
                        if (a !== n)
                            for (var o = 0, s = (r = n = a).length; o < s; ++o)
                                if (r[o].name === e) {
                                    (r = r.slice()).splice(o, 1);
                                    break;
                                }
                        i.tween = r;
                    };
                }
                function X(t, e, n) {
                    var r, i;
                    if ("function" != typeof n) throw new Error();
                    return function () {
                        var a = D(this, t),
                            o = a.tween;
                        if (o !== r) {
                            i = (r = o).slice();
                            for (var s = { name: e, value: n }, l = 0, c = i.length; l < c; ++l)
                                if (i[l].name === e) {
                                    i[l] = s;
                                    break;
                                }
                            l === c && i.push(s);
                        }
                        a.tween = i;
                    };
                }
                function G(t, e, n) {
                    var r = t._id;
                    return (
                        t.each(function () {
                            var t = D(this, r);
                            (t.value || (t.value = {}))[e] = n.apply(this, arguments);
                        }),
                        function (t) {
                            return I(t, r).value[e];
                        }
                    );
                }
                var K = n(650),
                    Q = n(765),
                    J = n(843);
                function tt(t, e) {
                    var n;
                    return (
                        "number" == typeof e
                            ? q.Z
                            : e instanceof K.ZP
                              ? Q.ZP
                              : (n = (0, K.ZP)(e))
                                ? ((e = n), Q.ZP)
                                : J.Z
                    )(t, e);
                }
                function et(t) {
                    return function () {
                        this.removeAttribute(t);
                    };
                }
                function nt(t) {
                    return function () {
                        this.removeAttributeNS(t.space, t.local);
                    };
                }
                function rt(t, e, n) {
                    var r,
                        i,
                        a = n + "";
                    return function () {
                        var o = this.getAttribute(t);
                        return o === a ? null : o === r ? i : (i = e((r = o), n));
                    };
                }
                function it(t, e, n) {
                    var r,
                        i,
                        a = n + "";
                    return function () {
                        var o = this.getAttributeNS(t.space, t.local);
                        return o === a ? null : o === r ? i : (i = e((r = o), n));
                    };
                }
                function at(t, e, n) {
                    var r, i, a;
                    return function () {
                        var o,
                            s,
                            l = n(this);
                        if (null != l)
                            return (o = this.getAttribute(t)) === (s = l + "")
                                ? null
                                : o === r && s === i
                                  ? a
                                  : ((i = s), (a = e((r = o), l)));
                        this.removeAttribute(t);
                    };
                }
                function ot(t, e, n) {
                    var r, i, a;
                    return function () {
                        var o,
                            s,
                            l = n(this);
                        if (null != l)
                            return (o = this.getAttributeNS(t.space, t.local)) === (s = l + "")
                                ? null
                                : o === r && s === i
                                  ? a
                                  : ((i = s), (a = e((r = o), l)));
                        this.removeAttributeNS(t.space, t.local);
                    };
                }
                function st(t, e) {
                    var n, r;
                    function i() {
                        var i = e.apply(this, arguments);
                        return (
                            i !== r &&
                                (n =
                                    (r = i) &&
                                    (function (t, e) {
                                        return function (n) {
                                            this.setAttributeNS(t.space, t.local, e.call(this, n));
                                        };
                                    })(t, i)),
                            n
                        );
                    }
                    return ((i._value = e), i);
                }
                function lt(t, e) {
                    var n, r;
                    function i() {
                        var i = e.apply(this, arguments);
                        return (
                            i !== r &&
                                (n =
                                    (r = i) &&
                                    (function (t, e) {
                                        return function (n) {
                                            this.setAttribute(t, e.call(this, n));
                                        };
                                    })(t, i)),
                            n
                        );
                    }
                    return ((i._value = e), i);
                }
                function ct(t, e) {
                    return function () {
                        P(this, t).delay = +e.apply(this, arguments);
                    };
                }
                function ut(t, e) {
                    return (
                        (e = +e),
                        function () {
                            P(this, t).delay = e;
                        }
                    );
                }
                function ht(t, e) {
                    return function () {
                        D(this, t).duration = +e.apply(this, arguments);
                    };
                }
                function dt(t, e) {
                    return (
                        (e = +e),
                        function () {
                            D(this, t).duration = e;
                        }
                    );
                }
                var gt = a.selection.prototype.constructor;
                function ft(t) {
                    return function () {
                        this.style.removeProperty(t);
                    };
                }
                var pt = 0;
                function mt(t, e, n, r) {
                    ((this._groups = t), (this._parents = e), (this._name = n), (this._id = r));
                }
                function vt(t) {
                    return (0, a.selection)().transition(t);
                }
                function yt() {
                    return ++pt;
                }
                var bt = a.selection.prototype;
                ((mt.prototype = vt.prototype =
                    {
                        constructor: mt,
                        select: function (t) {
                            var e = this._name,
                                n = this._id;
                            "function" != typeof t && (t = (0, a.selector)(t));
                            for (var r = this._groups, i = r.length, o = new Array(i), s = 0; s < i; ++s)
                                for (var l, c, u = r[s], h = u.length, d = (o[s] = new Array(h)), g = 0; g < h; ++g)
                                    (l = u[g]) &&
                                        (c = t.call(l, l.__data__, g, u)) &&
                                        ("__data__" in l && (c.__data__ = l.__data__),
                                        (d[g] = c),
                                        L(d[g], e, n, g, d, I(l, n)));
                            return new mt(o, this._parents, e, n);
                        },
                        selectAll: function (t) {
                            var e = this._name,
                                n = this._id;
                            "function" != typeof t && (t = (0, a.selectorAll)(t));
                            for (var r = this._groups, i = r.length, o = [], s = [], l = 0; l < i; ++l)
                                for (var c, u = r[l], h = u.length, d = 0; d < h; ++d)
                                    if ((c = u[d])) {
                                        for (
                                            var g, f = t.call(c, c.__data__, d, u), p = I(c, n), m = 0, v = f.length;
                                            m < v;
                                            ++m
                                        )
                                            (g = f[m]) && L(g, e, n, m, f, p);
                                        (o.push(f), s.push(c));
                                    }
                            return new mt(o, s, e, n);
                        },
                        filter: function (t) {
                            "function" != typeof t && (t = (0, a.matcher)(t));
                            for (var e = this._groups, n = e.length, r = new Array(n), i = 0; i < n; ++i)
                                for (var o, s = e[i], l = s.length, c = (r[i] = []), u = 0; u < l; ++u)
                                    (o = s[u]) && t.call(o, o.__data__, u, s) && c.push(o);
                            return new mt(r, this._parents, this._name, this._id);
                        },
                        merge: function (t) {
                            if (t._id !== this._id) throw new Error();
                            for (
                                var e = this._groups,
                                    n = t._groups,
                                    r = e.length,
                                    i = n.length,
                                    a = Math.min(r, i),
                                    o = new Array(r),
                                    s = 0;
                                s < a;
                                ++s
                            )
                                for (
                                    var l, c = e[s], u = n[s], h = c.length, d = (o[s] = new Array(h)), g = 0;
                                    g < h;
                                    ++g
                                )
                                    (l = c[g] || u[g]) && (d[g] = l);
                            for (; s < r; ++s) o[s] = e[s];
                            return new mt(o, this._parents, this._name, this._id);
                        },
                        selection: function () {
                            return new gt(this._groups, this._parents);
                        },
                        transition: function () {
                            for (
                                var t = this._name, e = this._id, n = yt(), r = this._groups, i = r.length, a = 0;
                                a < i;
                                ++a
                            )
                                for (var o, s = r[a], l = s.length, c = 0; c < l; ++c)
                                    if ((o = s[c])) {
                                        var u = I(o, e);
                                        L(o, t, n, c, s, {
                                            time: u.time + u.delay + u.duration,
                                            delay: 0,
                                            duration: u.duration,
                                            ease: u.ease,
                                        });
                                    }
                            return new mt(r, this._parents, t, n);
                        },
                        call: bt.call,
                        nodes: bt.nodes,
                        node: bt.node,
                        size: bt.size,
                        empty: bt.empty,
                        each: bt.each,
                        on: function (t, e) {
                            var n = this._id;
                            return arguments.length < 2
                                ? I(this.node(), n).on.on(t)
                                : this.each(
                                      (function (t, e, n) {
                                          var r,
                                              i,
                                              a = (function (t) {
                                                  return (t + "")
                                                      .trim()
                                                      .split(/^|\s+/)
                                                      .every(function (t) {
                                                          var e = t.indexOf(".");
                                                          return (e >= 0 && (t = t.slice(0, e)), !t || "start" === t);
                                                      });
                                              })(e)
                                                  ? P
                                                  : D;
                                          return function () {
                                              var o = a(this, t),
                                                  s = o.on;
                                              (s !== r && (i = (r = s).copy()).on(e, n), (o.on = i));
                                          };
                                      })(n, t, e),
                                  );
                        },
                        attr: function (t, e) {
                            var n = (0, a.namespace)(t),
                                r = "transform" === n ? Y : tt;
                            return this.attrTween(
                                t,
                                "function" == typeof e
                                    ? (n.local ? ot : at)(n, r, G(this, "attr." + t, e))
                                    : null == e
                                      ? (n.local ? nt : et)(n)
                                      : (n.local ? it : rt)(n, r, e),
                            );
                        },
                        attrTween: function (t, e) {
                            var n = "attr." + t;
                            if (arguments.length < 2) return (n = this.tween(n)) && n._value;
                            if (null == e) return this.tween(n, null);
                            if ("function" != typeof e) throw new Error();
                            var r = (0, a.namespace)(t);
                            return this.tween(n, (r.local ? st : lt)(r, e));
                        },
                        style: function (t, e, n) {
                            var r = "transform" == (t += "") ? W : tt;
                            return null == e
                                ? this.styleTween(
                                      t,
                                      (function (t, e) {
                                          var n, r, i;
                                          return function () {
                                              var o = (0, a.style)(this, t),
                                                  s = (this.style.removeProperty(t), (0, a.style)(this, t));
                                              return o === s
                                                  ? null
                                                  : o === n && s === r
                                                    ? i
                                                    : (i = e((n = o), (r = s)));
                                          };
                                      })(t, r),
                                  ).on("end.style." + t, ft(t))
                                : "function" == typeof e
                                  ? this.styleTween(
                                        t,
                                        (function (t, e, n) {
                                            var r, i, o;
                                            return function () {
                                                var s = (0, a.style)(this, t),
                                                    l = n(this),
                                                    c = l + "";
                                                return (
                                                    null == l &&
                                                        (this.style.removeProperty(t), (c = l = (0, a.style)(this, t))),
                                                    s === c
                                                        ? null
                                                        : s === r && c === i
                                                          ? o
                                                          : ((i = c), (o = e((r = s), l)))
                                                );
                                            };
                                        })(t, r, G(this, "style." + t, e)),
                                    ).each(
                                        (function (t, e) {
                                            var n,
                                                r,
                                                i,
                                                a,
                                                o = "style." + e,
                                                s = "end." + o;
                                            return function () {
                                                var l = D(this, t),
                                                    c = l.on,
                                                    u = null == l.value[o] ? a || (a = ft(e)) : void 0;
                                                ((c === n && i === u) || (r = (n = c).copy()).on(s, (i = u)),
                                                    (l.on = r));
                                            };
                                        })(this._id, t),
                                    )
                                  : this.styleTween(
                                        t,
                                        (function (t, e, n) {
                                            var r,
                                                i,
                                                o = n + "";
                                            return function () {
                                                var s = (0, a.style)(this, t);
                                                return s === o ? null : s === r ? i : (i = e((r = s), n));
                                            };
                                        })(t, r, e),
                                        n,
                                    ).on("end.style." + t, null);
                        },
                        styleTween: function (t, e, n) {
                            var r = "style." + (t += "");
                            if (arguments.length < 2) return (r = this.tween(r)) && r._value;
                            if (null == e) return this.tween(r, null);
                            if ("function" != typeof e) throw new Error();
                            return this.tween(
                                r,
                                (function (t, e, n) {
                                    var r, i;
                                    function a() {
                                        var a = e.apply(this, arguments);
                                        return (
                                            a !== i &&
                                                (r =
                                                    (i = a) &&
                                                    (function (t, e, n) {
                                                        return function (r) {
                                                            this.style.setProperty(t, e.call(this, r), n);
                                                        };
                                                    })(t, a, n)),
                                            r
                                        );
                                    }
                                    return ((a._value = e), a);
                                })(t, e, null == n ? "" : n),
                            );
                        },
                        text: function (t) {
                            return this.tween(
                                "text",
                                "function" == typeof t
                                    ? (function (t) {
                                          return function () {
                                              var e = t(this);
                                              this.textContent = null == e ? "" : e;
                                          };
                                      })(G(this, "text", t))
                                    : (function (t) {
                                          return function () {
                                              this.textContent = t;
                                          };
                                      })(null == t ? "" : t + ""),
                            );
                        },
                        textTween: function (t) {
                            var e = "text";
                            if (arguments.length < 1) return (e = this.tween(e)) && e._value;
                            if (null == t) return this.tween(e, null);
                            if ("function" != typeof t) throw new Error();
                            return this.tween(
                                e,
                                (function (t) {
                                    var e, n;
                                    function r() {
                                        var r = t.apply(this, arguments);
                                        return (
                                            r !== n &&
                                                (e =
                                                    (n = r) &&
                                                    (function (t) {
                                                        return function (e) {
                                                            this.textContent = t.call(this, e);
                                                        };
                                                    })(r)),
                                            e
                                        );
                                    }
                                    return ((r._value = t), r);
                                })(t),
                            );
                        },
                        remove: function () {
                            return this.on(
                                "end.remove",
                                (function (t) {
                                    return function () {
                                        var e = this.parentNode;
                                        for (var n in this.__transition) if (+n !== t) return;
                                        e && e.removeChild(this);
                                    };
                                })(this._id),
                            );
                        },
                        tween: function (t, e) {
                            var n = this._id;
                            if (((t += ""), arguments.length < 2)) {
                                for (var r, i = I(this.node(), n).tween, a = 0, o = i.length; a < o; ++a)
                                    if ((r = i[a]).name === t) return r.value;
                                return null;
                            }
                            return this.each((null == e ? Z : X)(n, t, e));
                        },
                        delay: function (t) {
                            var e = this._id;
                            return arguments.length
                                ? this.each(("function" == typeof t ? ct : ut)(e, t))
                                : I(this.node(), e).delay;
                        },
                        duration: function (t) {
                            var e = this._id;
                            return arguments.length
                                ? this.each(("function" == typeof t ? ht : dt)(e, t))
                                : I(this.node(), e).duration;
                        },
                        ease: function (t) {
                            var e = this._id;
                            return arguments.length
                                ? this.each(
                                      (function (t, e) {
                                          if ("function" != typeof e) throw new Error();
                                          return function () {
                                              D(this, t).ease = e;
                                          };
                                      })(e, t),
                                  )
                                : I(this.node(), e).ease;
                        },
                        end: function () {
                            var t,
                                e,
                                n = this,
                                r = n._id,
                                i = n.size();
                            return new Promise(function (a, o) {
                                var s = { value: o },
                                    l = {
                                        value: function () {
                                            0 == --i && a();
                                        },
                                    };
                                n.each(function () {
                                    var n = D(this, r),
                                        i = n.on;
                                    (i !== t &&
                                        ((e = (t = i).copy())._.cancel.push(s), e._.interrupt.push(s), e._.end.push(l)),
                                        (n.on = e));
                                });
                            });
                        },
                    }),
                    (function t(e) {
                        function n(t) {
                            return Math.pow(t, e);
                        }
                        return ((e = +e), (n.exponent = t), n);
                    })(3),
                    (function t(e) {
                        function n(t) {
                            return 1 - Math.pow(1 - t, e);
                        }
                        return ((e = +e), (n.exponent = t), n);
                    })(3),
                    (function t(e) {
                        function n(t) {
                            return ((t *= 2) <= 1 ? Math.pow(t, e) : 2 - Math.pow(2 - t, e)) / 2;
                        }
                        return ((e = +e), (n.exponent = t), n);
                    })(3),
                    Math.PI);
                var xt = 1.70158,
                    wt =
                        ((function t(e) {
                            function n(t) {
                                return t * t * ((e + 1) * t - e);
                            }
                            return ((e = +e), (n.overshoot = t), n);
                        })(xt),
                        (function t(e) {
                            function n(t) {
                                return --t * t * ((e + 1) * t + e) + 1;
                            }
                            return ((e = +e), (n.overshoot = t), n);
                        })(xt),
                        (function t(e) {
                            function n(t) {
                                return (
                                    ((t *= 2) < 1 ? t * t * ((e + 1) * t - e) : (t -= 2) * t * ((e + 1) * t + e) + 2) /
                                    2
                                );
                            }
                            return ((e = +e), (n.overshoot = t), n);
                        })(xt),
                        2 * Math.PI),
                    $t =
                        ((function t(e, n) {
                            var r = Math.asin(1 / (e = Math.max(1, e))) * (n /= wt);
                            function i(t) {
                                return e * Math.pow(2, 10 * --t) * Math.sin((r - t) / n);
                            }
                            return (
                                (i.amplitude = function (e) {
                                    return t(e, n * wt);
                                }),
                                (i.period = function (n) {
                                    return t(e, n);
                                }),
                                i
                            );
                        })(1, 0.3),
                        (function t(e, n) {
                            var r = Math.asin(1 / (e = Math.max(1, e))) * (n /= wt);
                            function i(t) {
                                return 1 - e * Math.pow(2, -10 * (t = +t)) * Math.sin((t + r) / n);
                            }
                            return (
                                (i.amplitude = function (e) {
                                    return t(e, n * wt);
                                }),
                                (i.period = function (n) {
                                    return t(e, n);
                                }),
                                i
                            );
                        })(1, 0.3),
                        (function t(e, n) {
                            var r = Math.asin(1 / (e = Math.max(1, e))) * (n /= wt);
                            function i(t) {
                                return (
                                    ((t = 2 * t - 1) < 0
                                        ? e * Math.pow(2, 10 * t) * Math.sin((r - t) / n)
                                        : 2 - e * Math.pow(2, -10 * t) * Math.sin((r + t) / n)) / 2
                                );
                            }
                            return (
                                (i.amplitude = function (e) {
                                    return t(e, n * wt);
                                }),
                                (i.period = function (n) {
                                    return t(e, n);
                                }),
                                i
                            );
                        })(1, 0.3),
                        {
                            time: null,
                            delay: 0,
                            duration: 250,
                            ease: function (t) {
                                return ((t *= 2) <= 1 ? t * t * t : (t -= 2) * t * t + 2) / 2;
                            },
                        });
                function At(t, e) {
                    for (var n; !(n = t.__transition) || !(n = n[e]); )
                        if (!(t = t.parentNode)) return (($t.time = m()), $t);
                    return n;
                }
                ((a.selection.prototype.interrupt = function (t) {
                    return this.each(function () {
                        E(this, t);
                    });
                }),
                    (a.selection.prototype.transition = function (t) {
                        var e, n;
                        t instanceof mt
                            ? ((e = t._id), (t = t._name))
                            : ((e = yt()), ((n = $t).time = m()), (t = null == t ? null : t + ""));
                        for (var r = this._groups, i = r.length, a = 0; a < i; ++a)
                            for (var o, s = r[a], l = s.length, c = 0; c < l; ++c)
                                (o = s[c]) && L(o, t, e, c, s, n || At(o, e));
                        return new mt(r, this._parents, t, e);
                    }));
                var _t = [null];
                function Ct(t, e) {
                    var n,
                        r,
                        i = t.__transition;
                    if (i)
                        for (r in ((e = null == e ? null : e + ""), i))
                            if ((n = i[r]).state > M && n.name === e) return new mt([[t]], _t, e, +r);
                    return null;
                }
            },
            261: (t, e, n) => {
                const r = n(881),
                    i = n(19);
                var a = !("undefined" == typeof window || !window.document || !window.document.createElement),
                    o = { Version: "1.0.0" };
                ((o.Drawer = r),
                    (o.Parser = i),
                    (o.clean = function (t) {
                        return t.replace(/[^A-Za-z0-9@\.\+\-\?!\(\)\[\]\{\}/\\=#\$:\*]/g, "");
                    }),
                    (o.apply = function (t, e = "canvas[data-smiles]", n = "light", i = null) {
                        let a = new r(t),
                            s = document.querySelectorAll(e);
                        for (var l = 0; l < s.length; l++) {
                            let t = s[l];
                            o.parse(
                                t.getAttribute("data-smiles"),
                                function (e) {
                                    a.draw(e, t, n, !1);
                                },
                                function (t) {
                                    i && i(t);
                                },
                            );
                        }
                    }),
                    (o.parse = function (t, e, n) {
                        try {
                            e && e(i.parse(t));
                        } catch (t) {
                            n && n(t);
                        }
                    }),
                    a && (window.SmilesDrawer = o),
                    Array.prototype.fill ||
                        Object.defineProperty(Array.prototype, "fill", {
                            value: function (t) {
                                if (null == this) throw new TypeError("this is null or not defined");
                                for (
                                    var e = Object(this),
                                        n = e.length >>> 0,
                                        r = arguments[1] >> 0,
                                        i = r < 0 ? Math.max(n + r, 0) : Math.min(r, n),
                                        a = arguments[2],
                                        o = void 0 === a ? n : a >> 0,
                                        s = o < 0 ? Math.max(n + o, 0) : Math.min(o, n);
                                    i < s;
                                )
                                    ((e[i] = t), i++);
                                return e;
                            },
                        }),
                    (t.exports = o));
            },
            348: (t) => {
                class e {
                    static clone(t) {
                        let n = Array.isArray(t) ? Array() : {};
                        for (let r in t) {
                            let i = t[r];
                            "function" == typeof i.clone
                                ? (n[r] = i.clone())
                                : (n[r] = "object" == typeof i ? e.clone(i) : i);
                        }
                        return n;
                    }
                    static equals(t, e) {
                        if (t.length !== e.length) return !1;
                        let n = t.slice().sort(),
                            r = e.slice().sort();
                        for (var i = 0; i < n.length; i++) if (n[i] !== r[i]) return !1;
                        return !0;
                    }
                    static print(t) {
                        if (0 == t.length) return "";
                        let e = "(";
                        for (let n = 0; n < t.length; n++) e += t[n].id ? t[n].id + ", " : t[n] + ", ";
                        return ((e = e.substring(0, e.length - 2)), e + ")");
                    }
                    static each(t, e) {
                        for (let n = 0; n < t.length; n++) e(t[n]);
                    }
                    static get(t, e, n) {
                        for (let r = 0; r < t.length; r++) if (t[r][e] == n) return t[r];
                    }
                    static contains(t, e) {
                        if (e.property || e.func) {
                            if (e.func) {
                                for (let n = 0; n < t.length; n++) if (e.func(t[n])) return !0;
                            } else for (let n = 0; n < t.length; n++) if (t[n][e.property] == e.value) return !0;
                        } else for (let n = 0; n < t.length; n++) if (t[n] == e.value) return !0;
                        return !1;
                    }
                    static intersection(t, e) {
                        let n = new Array();
                        for (let r = 0; r < t.length; r++)
                            for (let i = 0; i < e.length; i++) t[r] === e[i] && n.push(t[r]);
                        return n;
                    }
                    static unique(t) {
                        let e = {};
                        return t.filter(function (t) {
                            return void 0 === e[t] && (e[t] = !0);
                        });
                    }
                    static count(t, e) {
                        let n = 0;
                        for (let r = 0; r < t.length; r++) t[r] === e && n++;
                        return n;
                    }
                    static toggle(t, e) {
                        let n = Array(),
                            r = !1;
                        for (let i = 0; i < t.length; i++) t[i] !== e ? n.push(t[i]) : (r = !0);
                        return (r || n.push(e), n);
                    }
                    static remove(t, e) {
                        let n = Array();
                        for (let r = 0; r < t.length; r++) t[r] !== e && n.push(t[r]);
                        return n;
                    }
                    static removeUnique(t, e) {
                        let n = t.indexOf(e);
                        return (n > -1 && t.splice(n, 1), t);
                    }
                    static removeAll(t, e) {
                        return t.filter(function (t) {
                            return -1 === e.indexOf(t);
                        });
                    }
                    static merge(t, e) {
                        let n = new Array(t.length + e.length);
                        for (let e = 0; e < t.length; e++) n[e] = t[e];
                        for (let r = 0; r < e.length; r++) n[t.length + r] = e[r];
                        return n;
                    }
                    static containsAll(t, e) {
                        let n = 0;
                        for (let r = 0; r < t.length; r++) for (let i = 0; i < e.length; i++) t[r] === e[i] && n++;
                        return n === e.length;
                    }
                    static sortByAtomicNumberDesc(t) {
                        let e = t.map(function (t, e) {
                            return { index: e, value: t.atomicNumber.split(".").map(Number) };
                        });
                        return (
                            e.sort(function (t, e) {
                                let n = Math.min(e.value.length, t.value.length),
                                    r = 0;
                                for (; r < n && e.value[r] === t.value[r]; ) r++;
                                return r === n ? e.value.length - t.value.length : e.value[r] - t.value[r];
                            }),
                            e.map(function (e) {
                                return t[e.index];
                            })
                        );
                    }
                    static deepCopy(t) {
                        let n = Array();
                        for (let r = 0; r < t.length; r++) {
                            let i = t[r];
                            n[r] = i instanceof Array ? e.deepCopy(i) : i;
                        }
                        return n;
                    }
                }
                t.exports = e;
            },
            427: (t, e, n) => {
                const r = n(348);
                (n(152), n(421));
                class i {
                    constructor(t, e = "-") {
                        ((this.element = 1 === t.length ? t.toUpperCase() : t),
                            (this.drawExplicit = !1),
                            (this.ringbonds = Array()),
                            (this.rings = Array()),
                            (this.bondType = e),
                            (this.branchBond = null),
                            (this.isBridge = !1),
                            (this.isBridgeNode = !1),
                            (this.originalRings = Array()),
                            (this.bridgedRing = null),
                            (this.anchoredRings = Array()),
                            (this.bracket = null),
                            (this.plane = 0),
                            (this.attachedPseudoElements = {}),
                            (this.hasAttachedPseudoElements = !1),
                            (this.isDrawn = !0),
                            (this.isConnectedToRing = !1),
                            (this.neighbouringElements = Array()),
                            (this.isPartOfAromaticRing = t !== this.element),
                            (this.bondCount = 0),
                            (this.chirality = ""),
                            (this.isStereoCenter = !1),
                            (this.priority = 0),
                            (this.mainChain = !1),
                            (this.hydrogenDirection = "down"),
                            (this.subtreeDepth = 1),
                            (this.hasHydrogen = !1));
                    }
                    addNeighbouringElement(t) {
                        this.neighbouringElements.push(t);
                    }
                    attachPseudoElement(t, e, n = 0, r = 0) {
                        (null === n && (n = 0), null === r && (r = 0));
                        let i = n + t + r;
                        (this.attachedPseudoElements[i]
                            ? (this.attachedPseudoElements[i].count += 1)
                            : (this.attachedPseudoElements[i] = {
                                  element: t,
                                  count: 1,
                                  hydrogenCount: n,
                                  previousElement: e,
                                  charge: r,
                              }),
                            (this.hasAttachedPseudoElements = !0));
                    }
                    getAttachedPseudoElements() {
                        let t = {},
                            e = this;
                        return (
                            Object.keys(this.attachedPseudoElements)
                                .sort()
                                .forEach(function (n) {
                                    t[n] = e.attachedPseudoElements[n];
                                }),
                            t
                        );
                    }
                    getAttachedPseudoElementsCount() {
                        return Object.keys(this.attachedPseudoElements).length;
                    }
                    isHeteroAtom() {
                        return "C" !== this.element && "H" !== this.element;
                    }
                    addAnchoredRing(t) {
                        r.contains(this.anchoredRings, { value: t }) || this.anchoredRings.push(t);
                    }
                    getRingbondCount() {
                        return this.ringbonds.length;
                    }
                    backupRings() {
                        this.originalRings = Array(this.rings.length);
                        for (let t = 0; t < this.rings.length; t++) this.originalRings[t] = this.rings[t];
                    }
                    restoreRings() {
                        this.rings = Array(this.originalRings.length);
                        for (let t = 0; t < this.originalRings.length; t++) this.rings[t] = this.originalRings[t];
                    }
                    haveCommonRingbond(t, e) {
                        for (let n = 0; n < t.ringbonds.length; n++)
                            for (let r = 0; r < e.ringbonds.length; r++)
                                if (t.ringbonds[n].id == e.ringbonds[r].id) return !0;
                        return !1;
                    }
                    neighbouringElementsEqual(t) {
                        if (t.length !== this.neighbouringElements.length) return !1;
                        (t.sort(), this.neighbouringElements.sort());
                        for (var e = 0; e < this.neighbouringElements.length; e++)
                            if (t[e] !== this.neighbouringElements[e]) return !1;
                        return !0;
                    }
                    getAtomicNumber() {
                        return i.atomicNumbers[this.element];
                    }
                    getMaxBonds() {
                        return i.maxBonds[this.element];
                    }
                    static get maxBonds() {
                        return { H: 1, C: 4, N: 3, O: 2, P: 3, S: 2, B: 3, F: 1, I: 1, Cl: 1, Br: 1 };
                    }
                    static get atomicNumbers() {
                        return {
                            H: 1,
                            He: 2,
                            Li: 3,
                            Be: 4,
                            B: 5,
                            b: 5,
                            C: 6,
                            c: 6,
                            N: 7,
                            n: 7,
                            O: 8,
                            o: 8,
                            F: 9,
                            Ne: 10,
                            Na: 11,
                            Mg: 12,
                            Al: 13,
                            Si: 14,
                            P: 15,
                            p: 15,
                            S: 16,
                            s: 16,
                            Cl: 17,
                            Ar: 18,
                            K: 19,
                            Ca: 20,
                            Sc: 21,
                            Ti: 22,
                            V: 23,
                            Cr: 24,
                            Mn: 25,
                            Fe: 26,
                            Co: 27,
                            Ni: 28,
                            Cu: 29,
                            Zn: 30,
                            Ga: 31,
                            Ge: 32,
                            As: 33,
                            Se: 34,
                            Br: 35,
                            Kr: 36,
                            Rb: 37,
                            Sr: 38,
                            Y: 39,
                            Zr: 40,
                            Nb: 41,
                            Mo: 42,
                            Tc: 43,
                            Ru: 44,
                            Rh: 45,
                            Pd: 46,
                            Ag: 47,
                            Cd: 48,
                            In: 49,
                            Sn: 50,
                            Sb: 51,
                            Te: 52,
                            I: 53,
                            Xe: 54,
                            Cs: 55,
                            Ba: 56,
                            La: 57,
                            Ce: 58,
                            Pr: 59,
                            Nd: 60,
                            Pm: 61,
                            Sm: 62,
                            Eu: 63,
                            Gd: 64,
                            Tb: 65,
                            Dy: 66,
                            Ho: 67,
                            Er: 68,
                            Tm: 69,
                            Yb: 70,
                            Lu: 71,
                            Hf: 72,
                            Ta: 73,
                            W: 74,
                            Re: 75,
                            Os: 76,
                            Ir: 77,
                            Pt: 78,
                            Au: 79,
                            Hg: 80,
                            Tl: 81,
                            Pb: 82,
                            Bi: 83,
                            Po: 84,
                            At: 85,
                            Rn: 86,
                            Fr: 87,
                            Ra: 88,
                            Ac: 89,
                            Th: 90,
                            Pa: 91,
                            U: 92,
                            Np: 93,
                            Pu: 94,
                            Am: 95,
                            Cm: 96,
                            Bk: 97,
                            Cf: 98,
                            Es: 99,
                            Fm: 100,
                            Md: 101,
                            No: 102,
                            Lr: 103,
                            Rf: 104,
                            Db: 105,
                            Sg: 106,
                            Bh: 107,
                            Hs: 108,
                            Mt: 109,
                            Ds: 110,
                            Rg: 111,
                            Cn: 112,
                            Uut: 113,
                            Uuq: 114,
                            Uup: 115,
                            Uuh: 116,
                            Uus: 117,
                            Uuo: 118,
                        };
                    }
                    static get mass() {
                        return {
                            H: 1,
                            He: 2,
                            Li: 3,
                            Be: 4,
                            B: 5,
                            b: 5,
                            C: 6,
                            c: 6,
                            N: 7,
                            n: 7,
                            O: 8,
                            o: 8,
                            F: 9,
                            Ne: 10,
                            Na: 11,
                            Mg: 12,
                            Al: 13,
                            Si: 14,
                            P: 15,
                            p: 15,
                            S: 16,
                            s: 16,
                            Cl: 17,
                            Ar: 18,
                            K: 19,
                            Ca: 20,
                            Sc: 21,
                            Ti: 22,
                            V: 23,
                            Cr: 24,
                            Mn: 25,
                            Fe: 26,
                            Co: 27,
                            Ni: 28,
                            Cu: 29,
                            Zn: 30,
                            Ga: 31,
                            Ge: 32,
                            As: 33,
                            Se: 34,
                            Br: 35,
                            Kr: 36,
                            Rb: 37,
                            Sr: 38,
                            Y: 39,
                            Zr: 40,
                            Nb: 41,
                            Mo: 42,
                            Tc: 43,
                            Ru: 44,
                            Rh: 45,
                            Pd: 46,
                            Ag: 47,
                            Cd: 48,
                            In: 49,
                            Sn: 50,
                            Sb: 51,
                            Te: 52,
                            I: 53,
                            Xe: 54,
                            Cs: 55,
                            Ba: 56,
                            La: 57,
                            Ce: 58,
                            Pr: 59,
                            Nd: 60,
                            Pm: 61,
                            Sm: 62,
                            Eu: 63,
                            Gd: 64,
                            Tb: 65,
                            Dy: 66,
                            Ho: 67,
                            Er: 68,
                            Tm: 69,
                            Yb: 70,
                            Lu: 71,
                            Hf: 72,
                            Ta: 73,
                            W: 74,
                            Re: 75,
                            Os: 76,
                            Ir: 77,
                            Pt: 78,
                            Au: 79,
                            Hg: 80,
                            Tl: 81,
                            Pb: 82,
                            Bi: 83,
                            Po: 84,
                            At: 85,
                            Rn: 86,
                            Fr: 87,
                            Ra: 88,
                            Ac: 89,
                            Th: 90,
                            Pa: 91,
                            U: 92,
                            Np: 93,
                            Pu: 94,
                            Am: 95,
                            Cm: 96,
                            Bk: 97,
                            Cf: 98,
                            Es: 99,
                            Fm: 100,
                            Md: 101,
                            No: 102,
                            Lr: 103,
                            Rf: 104,
                            Db: 105,
                            Sg: 106,
                            Bh: 107,
                            Hs: 108,
                            Mt: 109,
                            Ds: 110,
                            Rg: 111,
                            Cn: 112,
                            Uut: 113,
                            Uuq: 114,
                            Uup: 115,
                            Uuh: 116,
                            Uus: 117,
                            Uuo: 118,
                        };
                    }
                }
                t.exports = i;
            },
            841: (t, e, n) => {
                const r = n(474),
                    i = n(614);
                (n(929),
                    n(152),
                    n(421),
                    (t.exports = class {
                        constructor(t, e, n) {
                            ((this.canvas =
                                "string" == typeof t || t instanceof String ? document.getElementById(t) : t),
                                (this.ctx = this.canvas.getContext("2d")),
                                (this.colors = e),
                                (this.opts = n),
                                (this.drawingWidth = 0),
                                (this.drawingHeight = 0),
                                (this.offsetX = 0),
                                (this.offsetY = 0),
                                (this.fontLarge = this.opts.fontSizeLarge + "pt Helvetica, Arial, sans-serif"),
                                (this.fontSmall = this.opts.fontSizeSmall + "pt Helvetica, Arial, sans-serif"),
                                this.updateSize(this.opts.width, this.opts.height),
                                (this.ctx.font = this.fontLarge),
                                (this.hydrogenWidth = this.ctx.measureText("H").width),
                                (this.halfHydrogenWidth = this.hydrogenWidth / 2),
                                (this.halfBondThickness = this.opts.bondThickness / 2));
                        }
                        updateSize(t, e) {
                            ((this.devicePixelRatio = window.devicePixelRatio || 1),
                                (this.backingStoreRatio =
                                    this.ctx.webkitBackingStorePixelRatio ||
                                    this.ctx.mozBackingStorePixelRatio ||
                                    this.ctx.msBackingStorePixelRatio ||
                                    this.ctx.oBackingStorePixelRatio ||
                                    this.ctx.backingStorePixelRatio ||
                                    1),
                                (this.ratio = this.devicePixelRatio / this.backingStoreRatio),
                                1 !== this.ratio
                                    ? ((this.canvas.width = t * this.ratio),
                                      (this.canvas.height = e * this.ratio),
                                      (this.canvas.style.width = t + "px"),
                                      (this.canvas.style.height = e + "px"),
                                      this.ctx.setTransform(this.ratio, 0, 0, this.ratio, 0, 0))
                                    : ((this.canvas.width = t * this.ratio), (this.canvas.height = e * this.ratio)));
                        }
                        setTheme(t) {
                            this.colors = t;
                        }
                        scale(t) {
                            let e = -Number.MAX_VALUE,
                                n = -Number.MAX_VALUE,
                                r = Number.MAX_VALUE,
                                i = Number.MAX_VALUE;
                            for (var a = 0; a < t.length; a++) {
                                if (!t[a].value.isDrawn) continue;
                                let o = t[a].position;
                                (e < o.x && (e = o.x),
                                    n < o.y && (n = o.y),
                                    r > o.x && (r = o.x),
                                    i > o.y && (i = o.y));
                            }
                            var o = this.opts.padding;
                            ((e += o),
                                (n += o),
                                (r -= o),
                                (i -= o),
                                (this.drawingWidth = e - r),
                                (this.drawingHeight = n - i));
                            var s = this.canvas.offsetWidth / this.drawingWidth,
                                l = this.canvas.offsetHeight / this.drawingHeight,
                                c = s < l ? s : l;
                            (this.ctx.scale(c, c),
                                (this.offsetX = -r),
                                (this.offsetY = -i),
                                s < l
                                    ? (this.offsetY += this.canvas.offsetHeight / (2 * c) - this.drawingHeight / 2)
                                    : (this.offsetX += this.canvas.offsetWidth / (2 * c) - this.drawingWidth / 2));
                        }
                        reset() {
                            this.ctx.setTransform(1, 0, 0, 1, 0, 0);
                        }
                        getColor(t) {
                            return (t = t.toUpperCase()) in this.colors ? this.colors[t] : this.colors.C;
                        }
                        drawCircle(t, e, n, i, a = !0, o = !1, s = "") {
                            let l = this.ctx,
                                c = this.offsetX,
                                u = this.offsetY;
                            (l.save(),
                                (l.lineWidth = 1.5),
                                l.beginPath(),
                                l.arc(t + c, e + u, n, 0, r.twoPI, !0),
                                l.closePath(),
                                o
                                    ? (a ? ((l.fillStyle = "#f00"), l.fill()) : ((l.strokeStyle = "#f00"), l.stroke()),
                                      this.drawDebugText(t, e, s))
                                    : a
                                      ? ((l.fillStyle = i), l.fill())
                                      : ((l.strokeStyle = i), l.stroke()),
                                l.restore());
                        }
                        drawLine(t, e = !1, n = 1) {
                            let r = this.ctx,
                                i = this.offsetX,
                                a = this.offsetY,
                                o = t.clone().shorten(4),
                                s = o.getLeftVector().clone(),
                                l = o.getRightVector().clone();
                            ((s.x += i),
                                (s.y += a),
                                (l.x += i),
                                (l.y += a),
                                e ||
                                    (r.save(),
                                    (r.globalCompositeOperation = "destination-out"),
                                    r.beginPath(),
                                    r.moveTo(s.x, s.y),
                                    r.lineTo(l.x, l.y),
                                    (r.lineCap = "round"),
                                    (r.lineWidth = this.opts.bondThickness + 1.2),
                                    (r.strokeStyle = this.getColor("BACKGROUND")),
                                    r.stroke(),
                                    (r.globalCompositeOperation = "source-over"),
                                    r.restore()),
                                (s = t.getLeftVector().clone()),
                                (l = t.getRightVector().clone()),
                                (s.x += i),
                                (s.y += a),
                                (l.x += i),
                                (l.y += a),
                                r.save(),
                                r.beginPath(),
                                r.moveTo(s.x, s.y),
                                r.lineTo(l.x, l.y),
                                (r.lineCap = "round"),
                                (r.lineWidth = this.opts.bondThickness));
                            let c = this.ctx.createLinearGradient(s.x, s.y, l.x, l.y);
                            (c.addColorStop(0.4, this.getColor(t.getLeftElement()) || this.getColor("C")),
                                c.addColorStop(0.6, this.getColor(t.getRightElement()) || this.getColor("C")),
                                e && (r.setLineDash([1, 1.5]), (r.lineWidth = this.opts.bondThickness / 1.5)),
                                n < 1 && (r.globalAlpha = n),
                                (r.strokeStyle = c),
                                r.stroke(),
                                r.restore());
                        }
                        drawWedge(t, e = 1) {
                            if (isNaN(t.from.x) || isNaN(t.from.y) || isNaN(t.to.x) || isNaN(t.to.y)) return;
                            let n = this.ctx,
                                r = this.offsetX,
                                a = this.offsetY,
                                o = t.clone().shorten(5),
                                s = o.getLeftVector().clone(),
                                l = o.getRightVector().clone();
                            ((s.x += r),
                                (s.y += a),
                                (l.x += r),
                                (l.y += a),
                                (s = t.getLeftVector().clone()),
                                (l = t.getRightVector().clone()),
                                (s.x += r),
                                (s.y += a),
                                (l.x += r),
                                (l.y += a),
                                n.save());
                            let c = i.normals(s, l);
                            (c[0].normalize(), c[1].normalize());
                            let u = s,
                                h = l;
                            t.getRightChiral() && ((u = l), (h = s));
                            let d = i.add(u, i.multiplyScalar(c[0], this.halfBondThickness)),
                                g = i.add(h, i.multiplyScalar(c[0], 1.5 + this.halfBondThickness)),
                                f = i.add(h, i.multiplyScalar(c[1], 1.5 + this.halfBondThickness)),
                                p = i.add(u, i.multiplyScalar(c[1], this.halfBondThickness));
                            (n.beginPath(),
                                n.moveTo(d.x, d.y),
                                n.lineTo(g.x, g.y),
                                n.lineTo(f.x, f.y),
                                n.lineTo(p.x, p.y));
                            let m = this.ctx.createRadialGradient(l.x, l.y, this.opts.bondLength, l.x, l.y, 0);
                            (m.addColorStop(0.4, this.getColor(t.getLeftElement()) || this.getColor("C")),
                                m.addColorStop(0.6, this.getColor(t.getRightElement()) || this.getColor("C")),
                                (n.fillStyle = m),
                                n.fill(),
                                n.restore());
                        }
                        drawDashedWedge(t) {
                            if (isNaN(t.from.x) || isNaN(t.from.y) || isNaN(t.to.x) || isNaN(t.to.y)) return;
                            let e = this.ctx,
                                n = this.offsetX,
                                r = this.offsetY,
                                a = t.getLeftVector().clone(),
                                o = t.getRightVector().clone();
                            ((a.x += n), (a.y += r), (o.x += n), (o.y += r), e.save());
                            let s = i.normals(a, o);
                            (s[0].normalize(), s[1].normalize());
                            let l,
                                c,
                                u,
                                h,
                                d = t.getRightChiral(),
                                g = t.clone();
                            (d
                                ? ((l = o),
                                  (c = a),
                                  g.shortenRight(1),
                                  (u = g.getRightVector().clone()),
                                  (h = g.getLeftVector().clone()))
                                : ((l = a),
                                  (c = o),
                                  g.shortenLeft(1),
                                  (u = g.getLeftVector().clone()),
                                  (h = g.getRightVector().clone())),
                                (u.x += n),
                                (u.y += r),
                                (h.x += n),
                                (h.y += r));
                            let f = i.subtract(c, l).normalize();
                            ((e.strokeStyle = this.getColor("C")),
                                (e.lineCap = "round"),
                                (e.lineWidth = this.opts.bondThickness),
                                e.beginPath());
                            let p = t.getLength(),
                                m = 1.25 / (p / (3 * this.opts.bondThickness)),
                                v = !1;
                            for (var y = 0; y < 1; y += m) {
                                let n = i.multiplyScalar(f, y * p),
                                    r = i.add(l, n),
                                    a = 1.5 * y,
                                    o = i.multiplyScalar(s[0], a);
                                (!v &&
                                    y > 0.5 &&
                                    (e.stroke(),
                                    e.beginPath(),
                                    (e.strokeStyle = this.getColor(t.getRightElement()) || this.getColor("C")),
                                    (v = !0)),
                                    r.subtract(o),
                                    e.moveTo(r.x, r.y),
                                    r.add(i.multiplyScalar(o, 2)),
                                    e.lineTo(r.x, r.y));
                            }
                            (e.stroke(), e.restore());
                        }
                        drawDebugText(t, e, n) {
                            let r = this.ctx;
                            (r.save(),
                                (r.font = "5px Droid Sans, sans-serif"),
                                (r.textAlign = "start"),
                                (r.textBaseline = "top"),
                                (r.fillStyle = "#ff0000"),
                                r.fillText(n, t + this.offsetX, e + this.offsetY),
                                r.restore());
                        }
                        drawBall(t, e, n) {
                            let i = this.ctx;
                            (i.save(),
                                i.beginPath(),
                                i.arc(t + this.offsetX, e + this.offsetY, this.opts.bondLength / 4.5, 0, r.twoPI, !1),
                                (i.fillStyle = this.getColor(n)),
                                i.fill(),
                                i.restore());
                        }
                        drawPoint(t, e, n) {
                            let i = this.ctx,
                                a = this.offsetX,
                                o = this.offsetY;
                            (i.save(),
                                (i.globalCompositeOperation = "destination-out"),
                                i.beginPath(),
                                i.arc(t + a, e + o, 1.5, 0, r.twoPI, !0),
                                i.closePath(),
                                i.fill(),
                                (i.globalCompositeOperation = "source-over"),
                                i.beginPath(),
                                i.arc(t + this.offsetX, e + this.offsetY, 0.75, 0, r.twoPI, !1),
                                (i.fillStyle = this.getColor(n)),
                                i.fill(),
                                i.restore());
                        }
                        drawText(t, e, n, i, a, o, s, l, c = {}) {
                            let u = this.ctx,
                                h = this.offsetX,
                                d = this.offsetY;
                            (u.save(), (u.textAlign = "start"), (u.textBaseline = "alphabetic"));
                            let g = "",
                                f = 0;
                            s && ((g = this.getChargeText(s)), (u.font = this.fontSmall), (f = u.measureText(g).width));
                            let p = "0",
                                m = 0;
                            (l > 0 && ((p = l.toString()), (u.font = this.fontSmall), (m = u.measureText(p).width)),
                                1 === s &&
                                    "N" === n &&
                                    c.hasOwnProperty("0O") &&
                                    c.hasOwnProperty("0O-1") &&
                                    ((c = {
                                        "0O": {
                                            element: "O",
                                            count: 2,
                                            hydrogenCount: 0,
                                            previousElement: "C",
                                            charge: "",
                                        },
                                    }),
                                    (s = 0)),
                                (u.font = this.fontLarge),
                                (u.fillStyle = this.getColor("BACKGROUND")));
                            let v = u.measureText(n);
                            ((v.totalWidth = v.width + f), (v.height = parseInt(this.fontLarge, 10)));
                            let y = v.width > this.opts.fontSizeLarge ? v.width : this.opts.fontSizeLarge;
                            ((y /= 1.5),
                                (u.globalCompositeOperation = "destination-out"),
                                u.beginPath(),
                                u.arc(t + h, e + d, y, 0, r.twoPI, !0),
                                u.closePath(),
                                u.fill(),
                                (u.globalCompositeOperation = "source-over"));
                            let b = -v.width / 2,
                                x = -v.width / 2;
                            ((u.fillStyle = this.getColor(n)),
                                u.fillText(n, t + h + b, e + this.opts.halfFontSizeLarge + d),
                                (b += v.width),
                                s &&
                                    ((u.font = this.fontSmall),
                                    u.fillText(g, t + h + b, e - this.opts.fifthFontSizeSmall + d),
                                    (b += f)),
                                l > 0 &&
                                    ((u.font = this.fontSmall),
                                    u.fillText(p, t + h + x - m, e - this.opts.fifthFontSizeSmall + d),
                                    (x -= m)),
                                (u.font = this.fontLarge));
                            let w = 0,
                                $ = 0;
                            if (1 === i) {
                                let n = t + h,
                                    r = e + d + this.opts.halfFontSizeLarge;
                                ((w = this.hydrogenWidth),
                                    (x -= w),
                                    "left" === a
                                        ? (n += x)
                                        : "right" === a || ("up" === a && o) || ("down" === a && o)
                                          ? (n += b)
                                          : "up" !== a || o
                                            ? "down" !== a ||
                                              o ||
                                              ((r += this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                                              (n -= this.halfHydrogenWidth))
                                            : ((r -= this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                                              (n -= this.halfHydrogenWidth)),
                                    u.fillText("H", n, r),
                                    (b += w));
                            } else if (i > 1) {
                                let n = t + h,
                                    r = e + d + this.opts.halfFontSizeLarge;
                                ((w = this.hydrogenWidth),
                                    (u.font = this.fontSmall),
                                    ($ = u.measureText(i).width),
                                    (x -= w + $),
                                    "left" === a
                                        ? (n += x)
                                        : "right" === a || ("up" === a && o) || ("down" === a && o)
                                          ? (n += b)
                                          : "up" !== a || o
                                            ? "down" !== a ||
                                              o ||
                                              ((r += this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                                              (n -= this.halfHydrogenWidth))
                                            : ((r -= this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                                              (n -= this.halfHydrogenWidth)),
                                    (u.font = this.fontLarge),
                                    u.fillText("H", n, r),
                                    (u.font = this.fontSmall),
                                    u.fillText(i, n + this.halfHydrogenWidth + $, r + this.opts.fifthFontSizeSmall),
                                    (b += w + this.halfHydrogenWidth + $));
                            }
                            for (let n in c) {
                                if (!c.hasOwnProperty(n)) continue;
                                let r = 0,
                                    i = 0,
                                    o = c[n].element,
                                    s = c[n].count,
                                    l = c[n].hydrogenCount,
                                    g = c[n].charge;
                                ((u.font = this.fontLarge),
                                    s > 1 && l > 0 && ((r = u.measureText("(").width), (i = u.measureText(")").width)));
                                let f = u.measureText(o).width,
                                    p = 0,
                                    m = "",
                                    v = 0;
                                ((w = 0),
                                    l > 0 && (w = this.hydrogenWidth),
                                    (u.font = this.fontSmall),
                                    s > 1 && (p = u.measureText(s).width),
                                    0 !== g && ((m = this.getChargeText(g)), (v = u.measureText(m).width)),
                                    ($ = 0),
                                    l > 1 && ($ = u.measureText(l).width),
                                    (u.font = this.fontLarge));
                                let y = t + h,
                                    A = e + d + this.opts.halfFontSizeLarge;
                                ((u.fillStyle = this.getColor(o)),
                                    s > 0 && (x -= p),
                                    s > 1 &&
                                        l > 0 &&
                                        ("left" === a
                                            ? ((x -= i), u.fillText(")", y + x, A))
                                            : (u.fillText("(", y + b, A), (b += r))),
                                    "left" === a
                                        ? ((x -= f), u.fillText(o, y + x, A))
                                        : (u.fillText(o, y + b, A), (b += f)),
                                    l > 0 &&
                                        ("left" === a
                                            ? ((x -= w + $),
                                              u.fillText("H", y + x, A),
                                              l > 1 &&
                                                  ((u.font = this.fontSmall),
                                                  u.fillText(l, y + x + w, A + this.opts.fifthFontSizeSmall)))
                                            : (u.fillText("H", y + b, A),
                                              (b += w),
                                              l > 1 &&
                                                  ((u.font = this.fontSmall),
                                                  u.fillText(l, y + b, A + this.opts.fifthFontSizeSmall),
                                                  (b += $)))),
                                    (u.font = this.fontLarge),
                                    s > 1 &&
                                        l > 0 &&
                                        ("left" === a
                                            ? ((x -= r), u.fillText("(", y + x, A))
                                            : (u.fillText(")", y + b, A), (b += i))),
                                    (u.font = this.fontSmall),
                                    s > 1 &&
                                        ("left" === a
                                            ? u.fillText(s, y + x + r + i + w + $ + f, A + this.opts.fifthFontSizeSmall)
                                            : (u.fillText(s, y + b, A + this.opts.fifthFontSizeSmall), (b += p))),
                                    0 !== g &&
                                        ("left" === a
                                            ? u.fillText(
                                                  m,
                                                  y + x + r + i + w + $ + f,
                                                  e - this.opts.fifthFontSizeSmall + d,
                                              )
                                            : (u.fillText(m, y + b, e - this.opts.fifthFontSizeSmall + d), (b += v))));
                            }
                            u.restore();
                        }
                        getChargeText(t) {
                            return 1 === t ? "+" : 2 === t ? "2+" : -1 === t ? "-" : -2 === t ? "2-" : "";
                        }
                        drawDebugPoint(t, e, n = "", r = "#f00") {
                            this.drawCircle(t, e, 2, r, !0, !0, n);
                        }
                        drawAromaticityRing(t) {
                            let e = this.ctx,
                                n = r.apothemFromSideLength(this.opts.bondLength, t.getSize());
                            (e.save(),
                                (e.strokeStyle = this.getColor("C")),
                                (e.lineWidth = this.opts.bondThickness),
                                e.beginPath(),
                                e.arc(
                                    t.center.x + this.offsetX,
                                    t.center.y + this.offsetY,
                                    n - this.opts.bondSpacing,
                                    0,
                                    2 * Math.PI,
                                    !0,
                                ),
                                e.closePath(),
                                e.stroke(),
                                e.restore());
                        }
                        clear() {
                            this.ctx.clearRect(0, 0, this.canvas.offsetWidth, this.canvas.offsetHeight);
                        }
                    }));
            },
            881: (t, e, n) => {
                const r = n(474),
                    i = n(348),
                    a = n(614),
                    o = n(929),
                    s = (n(152), n(826)),
                    l = n(427),
                    c = n(421),
                    u = n(333),
                    h = n(841),
                    d = n(707),
                    g = n(688);
                t.exports = class {
                    constructor(t) {
                        ((this.graph = null),
                            (this.doubleBondConfigCount = 0),
                            (this.doubleBondConfig = null),
                            (this.ringIdCounter = 0),
                            (this.ringConnectionIdCounter = 0),
                            (this.canvasWrapper = null),
                            (this.totalOverlapScore = 0),
                            (this.defaultOptions = {
                                width: 500,
                                height: 500,
                                bondThickness: 0.6,
                                bondLength: 15,
                                shortBondLength: 0.85,
                                bondSpacing: 0.18 * 15,
                                atomVisualization: "default",
                                isomeric: !0,
                                debug: !1,
                                terminalCarbons: !1,
                                explicitHydrogens: !1,
                                overlapSensitivity: 0.42,
                                overlapResolutionIterations: 1,
                                compactDrawing: !0,
                                fontSizeLarge: 5,
                                fontSizeSmall: 3,
                                padding: 20,
                                themes: {
                                    dark: {
                                        C: "#fff",
                                        O: "#e74c3c",
                                        N: "#3498db",
                                        F: "#27ae60",
                                        CL: "#16a085",
                                        BR: "#d35400",
                                        I: "#8e44ad",
                                        P: "#d35400",
                                        S: "#f1c40f",
                                        B: "#e67e22",
                                        SI: "#e67e22",
                                        H: "#fff",
                                        BACKGROUND: "#141414",
                                    },
                                    light: {
                                        C: "#222",
                                        O: "#e74c3c",
                                        N: "#3498db",
                                        F: "#27ae60",
                                        CL: "#16a085",
                                        BR: "#d35400",
                                        I: "#8e44ad",
                                        P: "#d35400",
                                        S: "#f1c40f",
                                        B: "#e67e22",
                                        SI: "#e67e22",
                                        H: "#222",
                                        BACKGROUND: "#fff",
                                    },
                                },
                            }),
                            (this.opts = this.extend(!0, this.defaultOptions, t)),
                            (this.opts.halfBondSpacing = this.opts.bondSpacing / 2),
                            (this.opts.bondLengthSq = this.opts.bondLength * this.opts.bondLength),
                            (this.opts.halfFontSizeLarge = this.opts.fontSizeLarge / 2),
                            (this.opts.quarterFontSizeLarge = this.opts.fontSizeLarge / 4),
                            (this.opts.fifthFontSizeSmall = this.opts.fontSizeSmall / 5),
                            (this.theme = this.opts.themes.dark));
                    }
                    extend() {
                        let t = this,
                            e = {},
                            n = !1,
                            r = 0,
                            i = arguments.length;
                        "[object Boolean]" === Object.prototype.toString.call(arguments[0]) &&
                            ((n = arguments[0]), r++);
                        let a = function (r) {
                            for (var i in r)
                                Object.prototype.hasOwnProperty.call(r, i) &&
                                    (n && "[object Object]" === Object.prototype.toString.call(r[i])
                                        ? (e[i] = t.extend(!0, e[i], r[i]))
                                        : (e[i] = r[i]));
                        };
                        for (; r < i; r++) a(arguments[r]);
                        return e;
                    }
                    draw(t, e, n = "light", i = !1) {
                        if (
                            ((this.data = t),
                            (this.infoOnly = i),
                            this.infoOnly || (this.canvasWrapper = new h(e, this.opts.themes[n], this.opts)),
                            (this.ringIdCounter = 0),
                            (this.ringConnectionIdCounter = 0),
                            (this.graph = new d(t, this.opts.isomeric)),
                            (this.rings = Array()),
                            (this.ringConnections = Array()),
                            (this.originalRings = Array()),
                            (this.originalRingConnections = Array()),
                            (this.bridgedRing = !1),
                            (this.doubleBondConfigCount = null),
                            (this.doubleBondConfig = null),
                            this.initRings(),
                            this.initHydrogens(),
                            !this.infoOnly)
                        ) {
                            (this.position(), this.restoreRingInformation(), this.resolvePrimaryOverlaps());
                            let t = this.getOverlapScore();
                            this.totalOverlapScore = this.getOverlapScore().total;
                            for (var a = 0; a < this.opts.overlapResolutionIterations; a++)
                                for (var o = 0; o < this.graph.edges.length; o++) {
                                    let e = this.graph.edges[o];
                                    if (this.isEdgeRotatable(e)) {
                                        let n = this.graph.getTreeDepth(e.sourceId, e.targetId),
                                            i = this.graph.getTreeDepth(e.targetId, e.sourceId),
                                            a = e.targetId,
                                            o = e.sourceId;
                                        if (
                                            (n > i && ((a = e.sourceId), (o = e.targetId)),
                                            this.getSubtreeOverlapScore(o, a, t.vertexScores).value >
                                                this.opts.overlapSensitivity)
                                        ) {
                                            let e = this.graph.vertices[a],
                                                n = this.graph.vertices[o],
                                                i = n.getNeighbours(a);
                                            if (1 === i.length) {
                                                let t = this.graph.vertices[i[0]],
                                                    a = t.position.getRotateAwayFromAngle(
                                                        e.position,
                                                        n.position,
                                                        r.toRad(120),
                                                    );
                                                this.rotateSubtree(t.id, n.id, a, n.position);
                                                let o = this.getOverlapScore().total;
                                                o > this.totalOverlapScore
                                                    ? this.rotateSubtree(t.id, n.id, -a, n.position)
                                                    : (this.totalOverlapScore = o);
                                            } else if (2 === i.length) {
                                                if (0 !== n.value.rings.length && 0 !== e.value.rings.length) continue;
                                                let t = this.graph.vertices[i[0]],
                                                    a = this.graph.vertices[i[1]];
                                                if (1 === t.value.rings.length && 1 === a.value.rings.length) {
                                                    if (t.value.rings[0] !== a.value.rings[0]) continue;
                                                } else {
                                                    if (0 !== t.value.rings.length || 0 !== a.value.rings.length)
                                                        continue;
                                                    {
                                                        let i = t.position.getRotateAwayFromAngle(
                                                                e.position,
                                                                n.position,
                                                                r.toRad(120),
                                                            ),
                                                            o = a.position.getRotateAwayFromAngle(
                                                                e.position,
                                                                n.position,
                                                                r.toRad(120),
                                                            );
                                                        (this.rotateSubtree(t.id, n.id, i, n.position),
                                                            this.rotateSubtree(a.id, n.id, o, n.position));
                                                        let s = this.getOverlapScore().total;
                                                        s > this.totalOverlapScore
                                                            ? (this.rotateSubtree(t.id, n.id, -i, n.position),
                                                              this.rotateSubtree(a.id, n.id, -o, n.position))
                                                            : (this.totalOverlapScore = s);
                                                    }
                                                }
                                            }
                                            t = this.getOverlapScore();
                                        }
                                    }
                                }
                            (this.resolveSecondaryOverlaps(t.scores),
                                this.opts.isomeric && this.annotateStereochemistry(),
                                this.opts.compactDrawing &&
                                    "default" === this.opts.atomVisualization &&
                                    this.initPseudoElements(),
                                this.rotateDrawing(),
                                this.canvasWrapper.scale(this.graph.vertices),
                                this.drawEdges(this.opts.debug),
                                this.drawVertices(this.opts.debug),
                                this.canvasWrapper.reset());
                        }
                    }
                    edgeRingCount(t) {
                        let e = this.graph.edges[t],
                            n = this.graph.vertices[e.sourceId],
                            r = this.graph.vertices[e.targetId];
                        return Math.min(n.value.rings.length, r.value.rings.length);
                    }
                    getBridgedRings() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isBridged && t.push(this.rings[e]);
                        return t;
                    }
                    getFusedRings() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isFused && t.push(this.rings[e]);
                        return t;
                    }
                    getSpiros() {
                        let t = Array();
                        for (var e = 0; e < this.rings.length; e++) this.rings[e].isSpiro && t.push(this.rings[e]);
                        return t;
                    }
                    printRingInfo() {
                        let t = "";
                        for (var e = 0; e < this.rings.length; e++) {
                            const n = this.rings[e];
                            ((t += n.id + ";"),
                                (t += n.members.length + ";"),
                                (t += n.neighbours.length + ";"),
                                (t += n.isSpiro ? "true;" : "false;"),
                                (t += n.isFused ? "true;" : "false;"),
                                (t += n.isBridged ? "true;" : "false;"),
                                (t += n.rings.length + ";"),
                                (t += "\n"));
                        }
                        return t;
                    }
                    rotateDrawing() {
                        let t = 0,
                            e = 0,
                            n = 0;
                        for (var r = 0; r < this.graph.vertices.length; r++) {
                            let a = this.graph.vertices[r];
                            if (a.value.isDrawn)
                                for (var i = r + 1; i < this.graph.vertices.length; i++) {
                                    let o = this.graph.vertices[i];
                                    if (!o.value.isDrawn) continue;
                                    let s = a.position.distanceSq(o.position);
                                    s > n && ((n = s), (t = r), (e = i));
                                }
                        }
                        let o = -a.subtract(this.graph.vertices[t].position, this.graph.vertices[e].position).angle();
                        if (!isNaN(o)) {
                            let t = o % 0.523599;
                            for (
                                t < 0.2617995 ? (o -= t) : (o += 0.523599 - t), r = 0;
                                r < this.graph.vertices.length;
                                r++
                            )
                                r !== e &&
                                    this.graph.vertices[r].position.rotateAround(o, this.graph.vertices[e].position);
                            for (r = 0; r < this.rings.length; r++)
                                this.rings[r].center.rotateAround(o, this.graph.vertices[e].position);
                        }
                    }
                    getTotalOverlapScore() {
                        return this.totalOverlapScore;
                    }
                    getRingCount() {
                        return this.rings.length;
                    }
                    hasBridgedRing() {
                        return this.bridgedRing;
                    }
                    getHeavyAtomCount() {
                        let t = 0;
                        for (var e = 0; e < this.graph.vertices.length; e++)
                            "H" !== this.graph.vertices[e].value.element && t++;
                        return t;
                    }
                    getMolecularFormula() {
                        let t = "",
                            e = new Map();
                        for (var n = 0; n < this.graph.vertices.length; n++) {
                            let t = this.graph.vertices[n].value;
                            if (
                                (e.has(t.element) ? e.set(t.element, e.get(t.element) + 1) : e.set(t.element, 1),
                                t.bracket &&
                                    !t.bracket.chirality &&
                                    (e.has("H")
                                        ? e.set("H", e.get("H") + t.bracket.hcount)
                                        : e.set("H", t.bracket.hcount)),
                                !t.bracket)
                            ) {
                                let n = l.maxBonds[t.element] - t.bondCount;
                                (t.isPartOfAromaticRing && n--,
                                    e.has("H") ? e.set("H", e.get("H") + n) : e.set("H", n));
                            }
                        }
                        if (e.has("C")) {
                            let n = e.get("C");
                            ((t += "C" + (n > 1 ? n : "")), e.delete("C"));
                        }
                        if (e.has("H")) {
                            let n = e.get("H");
                            ((t += "H" + (n > 1 ? n : "")), e.delete("H"));
                        }
                        return (
                            Object.keys(l.atomicNumbers)
                                .sort()
                                .map((n) => {
                                    if (e.has(n)) {
                                        let r = e.get(n);
                                        t += n + (r > 1 ? r : "");
                                    }
                                }),
                            t
                        );
                    }
                    getRingbondType(t, e) {
                        if (t.value.getRingbondCount() < 1 || e.value.getRingbondCount() < 1) return null;
                        for (var n = 0; n < t.value.ringbonds.length; n++)
                            for (var r = 0; r < e.value.ringbonds.length; r++)
                                if (t.value.ringbonds[n].id === e.value.ringbonds[r].id)
                                    return "-" === t.value.ringbonds[n].bondType
                                        ? e.value.ringbonds[r].bond
                                        : t.value.ringbonds[n].bond;
                        return null;
                    }
                    initRings() {
                        let t = new Map();
                        for (var e = this.graph.vertices.length - 1; e >= 0; e--) {
                            let r = this.graph.vertices[e];
                            if (0 !== r.value.ringbonds.length)
                                for (var n = 0; n < r.value.ringbonds.length; n++) {
                                    let e = r.value.ringbonds[n].id,
                                        i = r.value.ringbonds[n].bond;
                                    if (t.has(e)) {
                                        let a = r.id,
                                            o = t.get(e)[0],
                                            l = t.get(e)[1],
                                            c = new s(a, o, 1);
                                        c.setBondType(l || i || "-");
                                        let u = this.graph.addEdge(c),
                                            h = this.graph.vertices[o];
                                        (r.addRingbondChild(o, n),
                                            r.value.addNeighbouringElement(h.value.element),
                                            h.addRingbondChild(a, n),
                                            h.value.addNeighbouringElement(r.value.element),
                                            r.edges.push(u),
                                            h.edges.push(u),
                                            t.delete(e));
                                    } else t.set(e, [r.id, i]);
                                }
                        }
                        let r = g.getRings(this.graph);
                        if (null !== r) {
                            for (e = 0; e < r.length; e++) {
                                let t = [...r[e]],
                                    i = this.addRing(new c(t));
                                for (n = 0; n < t.length; n++) this.graph.vertices[t[n]].value.rings.push(i);
                            }
                            for (e = 0; e < this.rings.length - 1; e++)
                                for (n = e + 1; n < this.rings.length; n++) {
                                    let t = this.rings[e],
                                        r = this.rings[n],
                                        i = new u(t, r);
                                    i.vertices.size > 0 && this.addRingConnection(i);
                                }
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e];
                                t.neighbours = u.getNeighbours(this.ringConnections, t.id);
                            }
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e];
                                this.graph.vertices[t.members[0]].value.addAnchoredRing(t.id);
                            }
                            for (this.backupRingInformation(); this.rings.length > 0; ) {
                                let t = -1;
                                for (e = 0; e < this.rings.length; e++) {
                                    let n = this.rings[e];
                                    this.isPartOfBridgedRing(n.id) && !n.isBridged && (t = n.id);
                                }
                                if (-1 === t) break;
                                let n = this.getRing(t),
                                    r = this.getBridgedRingRings(n.id);
                                for (
                                    this.bridgedRing = !0, this.createBridgedRing(r, n.members[0]), e = 0;
                                    e < r.length;
                                    e++
                                )
                                    this.removeRing(r[e]);
                            }
                        }
                    }
                    initHydrogens() {
                        if (!this.opts.explicitHydrogens)
                            for (var t = 0; t < this.graph.vertices.length; t++) {
                                let e = this.graph.vertices[t];
                                if ("H" !== e.value.element) continue;
                                let n = this.graph.vertices[e.neighbours[0]];
                                ((n.value.hasHydrogen = !0),
                                    (!n.value.isStereoCenter ||
                                        (n.value.rings.length < 2 && !n.value.bridgedRing) ||
                                        (n.value.bridgedRing && n.value.originalRings.length < 2)) &&
                                        (e.value.isDrawn = !1));
                            }
                    }
                    getBridgedRingRings(t) {
                        let e = Array(),
                            n = this,
                            r = function (t) {
                                let i = n.getRing(t);
                                e.push(t);
                                for (var a = 0; a < i.neighbours.length; a++) {
                                    let o = i.neighbours[a];
                                    -1 === e.indexOf(o) &&
                                        o !== t &&
                                        u.isBridge(n.ringConnections, n.graph.vertices, t, o) &&
                                        r(o);
                                }
                            };
                        return (r(t), i.unique(e));
                    }
                    isPartOfBridgedRing(t) {
                        for (var e = 0; e < this.ringConnections.length; e++)
                            if (
                                this.ringConnections[e].containsRing(t) &&
                                this.ringConnections[e].isBridge(this.graph.vertices)
                            )
                                return !0;
                        return !1;
                    }
                    createBridgedRing(t, e) {
                        let n = new Set(),
                            r = new Set(),
                            a = new Set();
                        for (var o = 0; o < t.length; o++) {
                            let e = this.getRing(t[o]);
                            e.isPartOfBridged = !0;
                            for (var s = 0; s < e.members.length; s++) r.add(e.members[s]);
                            for (s = 0; s < e.neighbours.length; s++) {
                                let n = e.neighbours[s];
                                -1 === t.indexOf(n) && a.add(e.neighbours[s]);
                            }
                        }
                        let l = new Set();
                        for (let e of r) {
                            let r = this.graph.vertices[e],
                                a = i.intersection(t, r.value.rings);
                            1 === r.value.rings.length || 1 === a.length ? n.add(r.id) : l.add(r.id);
                        }
                        Array();
                        let u = Array();
                        for (let t of l) {
                            let e = this.graph.vertices[t],
                                r = !1;
                            for (let t = 0; t < e.edges.length; t++) 1 === this.edgeRingCount(e.edges[t]) && (r = !0);
                            r ? ((e.value.isBridgeNode = !0), n.add(e.id)) : ((e.value.isBridge = !0), n.add(e.id));
                        }
                        let h = new c([...n]);
                        for (h.isBridged = !0, h.neighbours = [...a], o = 0; o < t.length; o++)
                            h.rings.push(this.getRing(t[o]).clone());
                        for (this.addRing(h), o = 0; o < h.members.length; o++)
                            this.graph.vertices[h.members[o]].value.bridgedRing = h.id;
                        for (o = 0; o < u.length; o++) this.graph.vertices[u[o]].value.rings = Array();
                        for (let e of n) {
                            let n = this.graph.vertices[e];
                            ((n.value.rings = i.removeAll(n.value.rings, t)), n.value.rings.push(h.id));
                        }
                        for (o = 0; o < t.length; o++)
                            for (s = o + 1; s < t.length; s++) this.removeRingConnectionsBetween(t[o], t[s]);
                        for (let e of a) {
                            let n = this.getRingConnections(e, t);
                            for (s = 0; s < n.length; s++) this.getRingConnection(n[s]).updateOther(h.id, e);
                            this.getRing(e).neighbours.push(h.id);
                        }
                        return h;
                    }
                    areVerticesInSameRing(t, e) {
                        for (var n = 0; n < t.value.rings.length; n++)
                            for (var r = 0; r < e.value.rings.length; r++)
                                if (t.value.rings[n] === e.value.rings[r]) return !0;
                        return !1;
                    }
                    getCommonRings(t, e) {
                        let n = Array();
                        for (var r = 0; r < t.value.rings.length; r++)
                            for (var i = 0; i < e.value.rings.length; i++)
                                t.value.rings[r] == e.value.rings[i] && n.push(t.value.rings[r]);
                        return n;
                    }
                    getLargestOrAromaticCommonRing(t, e) {
                        let n = this.getCommonRings(t, e),
                            r = 0,
                            i = null;
                        for (var a = 0; a < n.length; a++) {
                            let t = this.getRing(n[a]),
                                e = t.getSize();
                            if (t.isBenzeneLike(this.graph.vertices)) return t;
                            e > r && ((r = e), (i = t));
                        }
                        return i;
                    }
                    getVerticesAt(t, e, n) {
                        let r = Array();
                        for (var i = 0; i < this.graph.vertices.length; i++) {
                            let a = this.graph.vertices[i];
                            a.id !== n && a.positioned && t.distanceSq(a.position) <= e * e && r.push(a.id);
                        }
                        return r;
                    }
                    getClosestVertex(t) {
                        let e = 99999,
                            n = null;
                        for (var r = 0; r < this.graph.vertices.length; r++) {
                            let i = this.graph.vertices[r];
                            if (i.id === t.id) continue;
                            let a = t.position.distanceSq(i.position);
                            a < e && ((e = a), (n = i));
                        }
                        return n;
                    }
                    addRing(t) {
                        return ((t.id = this.ringIdCounter++), this.rings.push(t), t.id);
                    }
                    removeRing(t) {
                        ((this.rings = this.rings.filter(function (e) {
                            return e.id !== t;
                        })),
                            (this.ringConnections = this.ringConnections.filter(function (e) {
                                return e.firstRingId !== t && e.secondRingId !== t;
                            })));
                        for (var e = 0; e < this.rings.length; e++) {
                            let n = this.rings[e];
                            n.neighbours = n.neighbours.filter(function (e) {
                                return e !== t;
                            });
                        }
                    }
                    getRing(t) {
                        for (var e = 0; e < this.rings.length; e++) if (this.rings[e].id == t) return this.rings[e];
                    }
                    addRingConnection(t) {
                        return ((t.id = this.ringConnectionIdCounter++), this.ringConnections.push(t), t.id);
                    }
                    removeRingConnection(t) {
                        this.ringConnections = this.ringConnections.filter(function (e) {
                            return e.id !== t;
                        });
                    }
                    removeRingConnectionsBetween(t, e) {
                        let n = Array();
                        for (var r = 0; r < this.ringConnections.length; r++) {
                            let i = this.ringConnections[r];
                            ((i.firstRingId === t && i.secondRingId === e) ||
                                (i.firstRingId === e && i.secondRingId === t)) &&
                                n.push(i.id);
                        }
                        for (r = 0; r < n.length; r++) this.removeRingConnection(n[r]);
                    }
                    getRingConnection(t) {
                        for (var e = 0; e < this.ringConnections.length; e++)
                            if (this.ringConnections[e].id == t) return this.ringConnections[e];
                    }
                    getRingConnections(t, e) {
                        let n = Array();
                        for (var r = 0; r < this.ringConnections.length; r++) {
                            let a = this.ringConnections[r];
                            for (var i = 0; i < e.length; i++) {
                                let r = e[i];
                                ((a.firstRingId === t && a.secondRingId === r) ||
                                    (a.firstRingId === r && a.secondRingId === t)) &&
                                    n.push(a.id);
                            }
                        }
                        return n;
                    }
                    getOverlapScore() {
                        let t = 0,
                            e = new Float32Array(this.graph.vertices.length);
                        for (var n = 0; n < this.graph.vertices.length; n++) e[n] = 0;
                        for (n = 0; n < this.graph.vertices.length; n++)
                            for (var r = this.graph.vertices.length; --r > n; ) {
                                let i = this.graph.vertices[n],
                                    o = this.graph.vertices[r];
                                if (!i.value.isDrawn || !o.value.isDrawn) continue;
                                let s = a.subtract(i.position, o.position).lengthSq();
                                if (s < this.opts.bondLengthSq) {
                                    let i = (this.opts.bondLength - Math.sqrt(s)) / this.opts.bondLength;
                                    ((t += i), (e[n] += i), (e[r] += i));
                                }
                            }
                        let i = Array();
                        for (n = 0; n < this.graph.vertices.length; n++) i.push({ id: n, score: e[n] });
                        return (
                            i.sort(function (t, e) {
                                return e.score - t.score;
                            }),
                            { total: t, scores: i, vertexScores: e }
                        );
                    }
                    chooseSide(t, e, n) {
                        let r = t.getNeighbours(e.id),
                            a = e.getNeighbours(t.id),
                            o = r.length,
                            s = a.length,
                            l = i.merge(r, a),
                            c = [0, 0];
                        for (var u = 0; u < l.length; u++)
                            this.graph.vertices[l[u]].position.sameSideAs(t.position, e.position, n[0])
                                ? c[0]++
                                : c[1]++;
                        let h = [0, 0];
                        for (u = 0; u < this.graph.vertices.length; u++)
                            this.graph.vertices[u].position.sameSideAs(t.position, e.position, n[0]) ? h[0]++ : h[1]++;
                        return {
                            totalSideCount: h,
                            totalPosition: h[0] > h[1] ? 0 : 1,
                            sideCount: c,
                            position: c[0] > c[1] ? 0 : 1,
                            anCount: o,
                            bnCount: s,
                        };
                    }
                    setRingCenter(t) {
                        let e = t.getSize(),
                            n = new a(0, 0);
                        for (var r = 0; r < e; r++) n.add(this.graph.vertices[t.members[r]].position);
                        t.center = n.divide(e);
                    }
                    getSubringCenter(t, e) {
                        let n = e.value.originalRings,
                            r = t.center,
                            i = Number.MAX_VALUE;
                        for (var a = 0; a < n.length; a++)
                            for (var o = 0; o < t.rings.length; o++)
                                n[a] === t.rings[o].id &&
                                    t.rings[o].getSize() < i &&
                                    ((r = t.rings[o].center), (i = t.rings[o].getSize()));
                        return r;
                    }
                    drawEdges(t) {
                        let e = this,
                            n = Array(this.graph.edges.length);
                        if (
                            (n.fill(!1),
                            this.graph.traverseBF(0, function (r) {
                                let i = e.graph.getEdges(r.id);
                                for (var a = 0; a < i.length; a++) {
                                    let r = i[a];
                                    n[r] || ((n[r] = !0), e.drawEdge(r, t));
                                }
                            }),
                            !this.bridgedRing)
                        )
                            for (var r = 0; r < this.rings.length; r++) {
                                let t = this.rings[r];
                                this.isRingAromatic(t) && this.canvasWrapper.drawAromaticityRing(t);
                            }
                    }
                    drawEdge(t, e) {
                        let n = this,
                            r = this.graph.edges[t],
                            s = this.graph.vertices[r.sourceId],
                            l = this.graph.vertices[r.targetId],
                            c = s.value.element,
                            u = l.value.element;
                        if (!((s.value.isDrawn && l.value.isDrawn) || "default" !== this.opts.atomVisualization))
                            return;
                        let h = s.position,
                            d = l.position,
                            g = this.getEdgeNormals(r),
                            f = i.clone(g);
                        if (
                            (f[0].multiplyScalar(10).add(h),
                            f[1].multiplyScalar(10).add(h),
                            "=" === r.bondType ||
                                "=" === this.getRingbondType(s, l) ||
                                (r.isPartOfAromaticRing && this.bridgedRing))
                        ) {
                            let t = this.areVerticesInSameRing(s, l),
                                e = this.chooseSide(s, l, f);
                            if (t) {
                                let t = this.getLargestOrAromaticCommonRing(s, l).center;
                                (g[0].multiplyScalar(n.opts.bondSpacing), g[1].multiplyScalar(n.opts.bondSpacing));
                                let e = null;
                                ((e = t.sameSideAs(s.position, l.position, a.add(h, g[0]))
                                    ? new o(a.add(h, g[0]), a.add(d, g[0]), c, u)
                                    : new o(a.add(h, g[1]), a.add(d, g[1]), c, u)),
                                    e.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength),
                                    r.isPartOfAromaticRing
                                        ? this.canvasWrapper.drawLine(e, !0)
                                        : this.canvasWrapper.drawLine(e),
                                    this.canvasWrapper.drawLine(new o(h, d, c, u)));
                            } else if (r.center || (s.isTerminal() && l.isTerminal())) {
                                (g[0].multiplyScalar(n.opts.halfBondSpacing),
                                    g[1].multiplyScalar(n.opts.halfBondSpacing));
                                let t = new o(a.add(h, g[0]), a.add(d, g[0]), c, u),
                                    e = new o(a.add(h, g[1]), a.add(d, g[1]), c, u);
                                (this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(e));
                            } else if ((0 == e.anCount && e.bnCount > 1) || (0 == e.bnCount && e.anCount > 1)) {
                                (g[0].multiplyScalar(n.opts.halfBondSpacing),
                                    g[1].multiplyScalar(n.opts.halfBondSpacing));
                                let t = new o(a.add(h, g[0]), a.add(d, g[0]), c, u),
                                    e = new o(a.add(h, g[1]), a.add(d, g[1]), c, u);
                                (this.canvasWrapper.drawLine(t), this.canvasWrapper.drawLine(e));
                            } else if (e.sideCount[0] > e.sideCount[1]) {
                                (g[0].multiplyScalar(n.opts.bondSpacing), g[1].multiplyScalar(n.opts.bondSpacing));
                                let t = new o(a.add(h, g[0]), a.add(d, g[0]), c, u);
                                (t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength),
                                    this.canvasWrapper.drawLine(t),
                                    this.canvasWrapper.drawLine(new o(h, d, c, u)));
                            } else if (e.sideCount[0] < e.sideCount[1]) {
                                (g[0].multiplyScalar(n.opts.bondSpacing), g[1].multiplyScalar(n.opts.bondSpacing));
                                let t = new o(a.add(h, g[1]), a.add(d, g[1]), c, u);
                                (t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength),
                                    this.canvasWrapper.drawLine(t),
                                    this.canvasWrapper.drawLine(new o(h, d, c, u)));
                            } else if (e.totalSideCount[0] > e.totalSideCount[1]) {
                                (g[0].multiplyScalar(n.opts.bondSpacing), g[1].multiplyScalar(n.opts.bondSpacing));
                                let t = new o(a.add(h, g[0]), a.add(d, g[0]), c, u);
                                (t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength),
                                    this.canvasWrapper.drawLine(t),
                                    this.canvasWrapper.drawLine(new o(h, d, c, u)));
                            } else if (e.totalSideCount[0] <= e.totalSideCount[1]) {
                                (g[0].multiplyScalar(n.opts.bondSpacing), g[1].multiplyScalar(n.opts.bondSpacing));
                                let t = new o(a.add(h, g[1]), a.add(d, g[1]), c, u);
                                (t.shorten(this.opts.bondLength - this.opts.shortBondLength * this.opts.bondLength),
                                    this.canvasWrapper.drawLine(t),
                                    this.canvasWrapper.drawLine(new o(h, d, c, u)));
                            }
                        } else if ("#" === r.bondType) {
                            (g[0].multiplyScalar(n.opts.bondSpacing / 1.5),
                                g[1].multiplyScalar(n.opts.bondSpacing / 1.5));
                            let t = new o(a.add(h, g[0]), a.add(d, g[0]), c, u),
                                e = new o(a.add(h, g[1]), a.add(d, g[1]), c, u);
                            (this.canvasWrapper.drawLine(t),
                                this.canvasWrapper.drawLine(e),
                                this.canvasWrapper.drawLine(new o(h, d, c, u)));
                        } else if ("." === r.bondType);
                        else {
                            let t = s.value.isStereoCenter,
                                e = l.value.isStereoCenter;
                            "up" === r.wedge
                                ? this.canvasWrapper.drawWedge(new o(h, d, c, u, t, e))
                                : "down" === r.wedge
                                  ? this.canvasWrapper.drawDashedWedge(new o(h, d, c, u, t, e))
                                  : this.canvasWrapper.drawLine(new o(h, d, c, u, t, e));
                        }
                        if (e) {
                            let e = a.midpoint(h, d);
                            this.canvasWrapper.drawDebugText(e.x, e.y, "e: " + t);
                        }
                    }
                    drawVertices(t) {
                        var e = this.graph.vertices.length;
                        for (e = 0; e < this.graph.vertices.length; e++) {
                            let n = this.graph.vertices[e],
                                r = n.value,
                                o = 0,
                                s = 0,
                                c = n.value.bondCount,
                                u = r.element,
                                h = l.maxBonds[u] - c,
                                d = n.getTextDirection(this.graph.vertices),
                                g =
                                    !(!this.opts.terminalCarbons && "C" === u && !r.hasAttachedPseudoElements) &&
                                    n.isTerminal(),
                                f = "C" === r.element;
                            if (
                                ("N" === r.element && r.isPartOfAromaticRing && (h = 0),
                                r.bracket && ((h = r.bracket.hcount), (o = r.bracket.charge), (s = r.bracket.isotope)),
                                "allballs" === this.opts.atomVisualization)
                            )
                                this.canvasWrapper.drawBall(n.position.x, n.position.y, u);
                            else if (
                                (r.isDrawn && (!f || r.drawExplicit || g || r.hasAttachedPseudoElements)) ||
                                1 === this.graph.vertices.length
                            )
                                "default" === this.opts.atomVisualization
                                    ? this.canvasWrapper.drawText(
                                          n.position.x,
                                          n.position.y,
                                          u,
                                          h,
                                          d,
                                          g,
                                          o,
                                          s,
                                          r.getAttachedPseudoElements(),
                                      )
                                    : "balls" === this.opts.atomVisualization &&
                                      this.canvasWrapper.drawBall(n.position.x, n.position.y, u);
                            else if (2 === n.getNeighbourCount() && 1 == n.forcePositioned) {
                                let t = this.graph.vertices[n.neighbours[0]].position,
                                    e = this.graph.vertices[n.neighbours[1]].position,
                                    r = a.threePointangle(n.position, t, e);
                                Math.abs(Math.PI - r) < 0.1 &&
                                    this.canvasWrapper.drawPoint(n.position.x, n.position.y, u);
                            }
                            if (t) {
                                let t = "v: " + n.id + " " + i.print(r.ringbonds);
                                this.canvasWrapper.drawDebugText(n.position.x, n.position.y, t);
                            }
                        }
                        if (this.opts.debug)
                            for (e = 0; e < this.rings.length; e++) {
                                let t = this.rings[e].center;
                                this.canvasWrapper.drawDebugPoint(t.x, t.y, "r: " + this.rings[e].id);
                            }
                    }
                    position() {
                        let t = null;
                        for (var e = 0; e < this.graph.vertices.length; e++)
                            if (null !== this.graph.vertices[e].value.bridgedRing) {
                                t = this.graph.vertices[e];
                                break;
                            }
                        for (e = 0; e < this.rings.length; e++)
                            this.rings[e].isBridged && (t = this.graph.vertices[this.rings[e].members[0]]);
                        (this.rings.length > 0 && null === t && (t = this.graph.vertices[this.rings[0].members[0]]),
                            null === t && (t = this.graph.vertices[0]),
                            this.createNextBond(t, null, 0));
                    }
                    backupRingInformation() {
                        ((this.originalRings = Array()), (this.originalRingConnections = Array()));
                        for (var t = 0; t < this.rings.length; t++) this.originalRings.push(this.rings[t]);
                        for (t = 0; t < this.ringConnections.length; t++)
                            this.originalRingConnections.push(this.ringConnections[t]);
                        for (t = 0; t < this.graph.vertices.length; t++) this.graph.vertices[t].value.backupRings();
                    }
                    restoreRingInformation() {
                        let t = this.getBridgedRings();
                        ((this.rings = Array()), (this.ringConnections = Array()));
                        for (var e = 0; e < t.length; e++) {
                            let r = t[e];
                            for (var n = 0; n < r.rings.length; n++) {
                                let t = r.rings[n];
                                this.originalRings[t.id].center = t.center;
                            }
                        }
                        for (e = 0; e < this.originalRings.length; e++) this.rings.push(this.originalRings[e]);
                        for (e = 0; e < this.originalRingConnections.length; e++)
                            this.ringConnections.push(this.originalRingConnections[e]);
                        for (e = 0; e < this.graph.vertices.length; e++) this.graph.vertices[e].value.restoreRings();
                    }
                    createRing(t, e = null, n = null, i = null) {
                        if (t.positioned) return;
                        e = e || new a(0, 0);
                        let o = t.getOrderedNeighbours(this.ringConnections),
                            s = n ? a.subtract(n.position, e).angle() : 0,
                            l = r.polyCircumradius(this.opts.bondLength, t.getSize()),
                            c = r.centralAngle(t.getSize());
                        t.centralAngle = c;
                        let h = s,
                            d = this,
                            g = n ? n.id : null;
                        if (
                            (-1 === t.members.indexOf(g) && (n && (n.positioned = !1), (g = t.members[0])), t.isBridged)
                        ) {
                            (this.graph.kkLayout(t.members.slice(), e, n.id, t, this.opts.bondLength),
                                (t.positioned = !0),
                                this.setRingCenter(t),
                                (e = t.center));
                            for (var f = 0; f < t.rings.length; f++) this.setRingCenter(t.rings[f]);
                        } else
                            t.eachMember(
                                this.graph.vertices,
                                function (n) {
                                    let r = d.graph.vertices[n];
                                    (r.positioned || r.setPosition(e.x + Math.cos(h) * l, e.y + Math.sin(h) * l),
                                        (h += c),
                                        (!t.isBridged || t.rings.length < 3) && ((r.angle = h), (r.positioned = !0)));
                                },
                                g,
                                i ? i.id : null,
                            );
                        for (t.positioned = !0, t.center = e, f = 0; f < o.length; f++) {
                            let n = this.getRing(o[f].neighbour);
                            if (n.positioned) continue;
                            let i = u.getVertices(this.ringConnections, t.id, n.id);
                            if (2 === i.length) {
                                ((t.isFused = !0), (n.isFused = !0));
                                let o = this.graph.vertices[i[0]],
                                    s = this.graph.vertices[i[1]],
                                    l = a.midpoint(o.position, s.position),
                                    c = a.normals(o.position, s.position);
                                (c[0].normalize(), c[1].normalize());
                                let u = r.polyCircumradius(this.opts.bondLength, n.getSize()),
                                    h = r.apothem(u, n.getSize());
                                (c[0].multiplyScalar(h).add(l), c[1].multiplyScalar(h).add(l));
                                let d = c[0];
                                a.subtract(e, c[1]).lengthSq() > a.subtract(e, c[0]).lengthSq() && (d = c[1]);
                                let g = a.subtract(o.position, d),
                                    f = a.subtract(s.position, d);
                                -1 === g.clockwise(f)
                                    ? n.positioned || this.createRing(n, d, o, s)
                                    : n.positioned || this.createRing(n, d, s, o);
                            } else if (1 === i.length) {
                                ((t.isSpiro = !0), (n.isSpiro = !0));
                                let o = this.graph.vertices[i[0]],
                                    s = a.subtract(e, o.position);
                                (s.invert(), s.normalize());
                                let l = r.polyCircumradius(this.opts.bondLength, n.getSize());
                                (s.multiplyScalar(l), s.add(o.position), n.positioned || this.createRing(n, s, o));
                            }
                        }
                        for (f = 0; f < t.members.length; f++) {
                            let e = this.graph.vertices[t.members[f]],
                                n = e.neighbours;
                            for (var p = 0; p < n.length; p++) {
                                let t = this.graph.vertices[n[p]];
                                t.positioned || ((t.value.isConnectedToRing = !0), this.createNextBond(t, e, 0));
                            }
                        }
                    }
                    rotateSubtree(t, e, n, r) {
                        let i = this;
                        this.graph.traverseTree(t, e, function (t) {
                            t.position.rotateAround(n, r);
                            for (var e = 0; e < t.value.anchoredRings.length; e++) {
                                let a = i.rings[t.value.anchoredRings[e]];
                                a && a.center.rotateAround(n, r);
                            }
                        });
                    }
                    getSubtreeOverlapScore(t, e, n) {
                        let r = this,
                            i = 0,
                            o = new a(0, 0),
                            s = 0;
                        return (
                            this.graph.traverseTree(t, e, function (t) {
                                if (!t.value.isDrawn) return;
                                let e = n[t.id];
                                e > r.opts.overlapSensitivity && ((i += e), s++);
                                let a = r.graph.vertices[t.id].position.clone();
                                (a.multiplyScalar(e), o.add(a));
                            }),
                            o.divide(i),
                            { value: i / s, center: o }
                        );
                    }
                    getCurrentCenterOfMass() {
                        let t = new a(0, 0),
                            e = 0;
                        for (var n = 0; n < this.graph.vertices.length; n++) {
                            let r = this.graph.vertices[n];
                            r.positioned && (t.add(r.position), e++);
                        }
                        return t.divide(e);
                    }
                    getCurrentCenterOfMassInNeigbourhood(t, e = 2 * this.opts.bondLength) {
                        let n = new a(0, 0),
                            r = 0,
                            i = e * e;
                        for (var o = 0; o < this.graph.vertices.length; o++) {
                            let e = this.graph.vertices[o];
                            e.positioned && t.distanceSq(e.position) < i && (n.add(e.position), r++);
                        }
                        return n.divide(r);
                    }
                    resolvePrimaryOverlaps() {
                        let t = Array(),
                            e = Array(this.graph.vertices.length);
                        for (var n = 0; n < this.rings.length; n++) {
                            let a = this.rings[n];
                            for (var r = 0; r < a.members.length; r++) {
                                let n = this.graph.vertices[a.members[r]];
                                if (e[n.id]) continue;
                                e[n.id] = !0;
                                let o = this.getNonRingNeighbours(n.id);
                                if (o.length > 1) {
                                    let e = Array();
                                    for (var i = 0; i < n.value.rings.length; i++) e.push(n.value.rings[i]);
                                    t.push({ common: n, rings: e, vertices: o });
                                } else if (1 === o.length && 2 === n.value.rings.length) {
                                    let e = Array();
                                    for (i = 0; i < n.value.rings.length; i++) e.push(n.value.rings[i]);
                                    t.push({ common: n, rings: e, vertices: o });
                                }
                            }
                        }
                        for (n = 0; n < t.length; n++) {
                            let e = t[n];
                            if (2 === e.vertices.length) {
                                let t = e.vertices[0],
                                    n = e.vertices[1];
                                if (!t.value.isDrawn || !n.value.isDrawn) continue;
                                let r = (2 * Math.PI - this.getRing(e.rings[0]).getAngle()) / 6;
                                (this.rotateSubtree(t.id, e.common.id, r, e.common.position),
                                    this.rotateSubtree(n.id, e.common.id, -r, e.common.position));
                                let i = this.getOverlapScore(),
                                    a = this.getSubtreeOverlapScore(t.id, e.common.id, i.vertexScores),
                                    o = this.getSubtreeOverlapScore(n.id, e.common.id, i.vertexScores),
                                    s = a.value + o.value;
                                (this.rotateSubtree(t.id, e.common.id, -2 * r, e.common.position),
                                    this.rotateSubtree(n.id, e.common.id, 2 * r, e.common.position),
                                    (i = this.getOverlapScore()),
                                    (a = this.getSubtreeOverlapScore(t.id, e.common.id, i.vertexScores)),
                                    (o = this.getSubtreeOverlapScore(n.id, e.common.id, i.vertexScores)),
                                    a.value + o.value > s &&
                                        (this.rotateSubtree(t.id, e.common.id, 2 * r, e.common.position),
                                        this.rotateSubtree(n.id, e.common.id, -2 * r, e.common.position)));
                            } else 1 === e.vertices.length && e.rings.length;
                        }
                    }
                    resolveSecondaryOverlaps(t) {
                        for (var e = 0; e < t.length; e++)
                            if (t[e].score > this.opts.overlapSensitivity) {
                                let n = this.graph.vertices[t[e].id];
                                if (n.isTerminal()) {
                                    let t = this.getClosestVertex(n);
                                    if (t) {
                                        let e = null;
                                        e = t.isTerminal()
                                            ? 0 === t.id
                                                ? this.graph.vertices[1].position
                                                : t.previousPosition
                                            : 0 === t.id
                                              ? this.graph.vertices[1].position
                                              : t.position;
                                        let i = 0 === n.id ? this.graph.vertices[1].position : n.previousPosition;
                                        n.position.rotateAwayFrom(e, i, r.toRad(20));
                                    }
                                }
                            }
                    }
                    getLastVertexWithAngle(t) {
                        let e = 0,
                            n = null;
                        for (; !e && t; ) ((n = this.graph.vertices[t]), (e = n.angle), (t = n.parentVertexId));
                        return n;
                    }
                    createNextBond(t, e = null, n = 0, o = !1, s = !1) {
                        if (t.positioned && !s) return;
                        let l = !1;
                        if (e) {
                            let n = this.graph.getEdge(t.id, e.id);
                            ("/" !== n.bondType && "\\" !== n.bondType) ||
                                ++this.doubleBondConfigCount % 2 != 1 ||
                                (null === this.doubleBondConfig &&
                                    ((this.doubleBondConfig = n.bondType),
                                    (l = !0),
                                    null === e.parentVertexId &&
                                        t.value.branchBond &&
                                        ("/" === this.doubleBondConfig
                                            ? (this.doubleBondConfig = "\\")
                                            : "\\" === this.doubleBondConfig && (this.doubleBondConfig = "/"))));
                        }
                        if (!s)
                            if (e)
                                if (e.value.rings.length > 0) {
                                    let n = e.neighbours,
                                        r = null,
                                        o = new a(0, 0);
                                    if (null === e.value.bridgedRing && e.value.rings.length > 1)
                                        for (var c = 0; c < n.length; c++) {
                                            let t = this.graph.vertices[n[c]];
                                            if (i.containsAll(t.value.rings, e.value.rings)) {
                                                r = t;
                                                break;
                                            }
                                        }
                                    if (null === r) {
                                        for (c = 0; c < n.length; c++) {
                                            let t = this.graph.vertices[n[c]];
                                            t.positioned &&
                                                this.areVerticesInSameRing(t, e) &&
                                                o.add(a.subtract(t.position, e.position));
                                        }
                                        o.invert().normalize().multiplyScalar(this.opts.bondLength).add(e.position);
                                    } else o = r.position.clone().rotateAround(Math.PI, e.position);
                                    ((t.previousPosition = e.position),
                                        t.setPositionFromVector(o),
                                        (t.positioned = !0));
                                } else {
                                    let r = new a(this.opts.bondLength, 0);
                                    (r.rotate(n),
                                        r.add(e.position),
                                        t.setPositionFromVector(r),
                                        (t.previousPosition = e.position),
                                        (t.positioned = !0));
                                }
                            else {
                                let e = new a(this.opts.bondLength, 0);
                                (e.rotate(r.toRad(-60)),
                                    (t.previousPosition = e),
                                    t.setPosition(this.opts.bondLength, 0),
                                    (t.angle = r.toRad(-60)),
                                    null === t.value.bridgedRing && (t.positioned = !0));
                            }
                        if (null !== t.value.bridgedRing) {
                            let e = this.getRing(t.value.bridgedRing);
                            if (!e.positioned) {
                                let n = a.subtract(t.previousPosition, t.position);
                                (n.invert(), n.normalize());
                                let i = r.polyCircumradius(this.opts.bondLength, e.members.length);
                                (n.multiplyScalar(i), n.add(t.position), this.createRing(e, n, t));
                            }
                        } else if (t.value.rings.length > 0) {
                            let e = this.getRing(t.value.rings[0]);
                            if (!e.positioned) {
                                let n = a.subtract(t.previousPosition, t.position);
                                (n.invert(), n.normalize());
                                let i = r.polyCircumradius(this.opts.bondLength, e.getSize());
                                (n.multiplyScalar(i), n.add(t.position), this.createRing(e, n, t));
                            }
                        } else {
                            t.value.isStereoCenter;
                            let n = t.getNeighbours(),
                                s = Array();
                            for (c = 0; c < n.length; c++) this.graph.vertices[n[c]].value.isDrawn && s.push(n[c]);
                            e && (s = i.remove(s, e.id));
                            let u = t.getAngle();
                            if (1 === s.length) {
                                let n = this.graph.vertices[s[0]];
                                if (
                                    "#" === t.value.bondType ||
                                    (e && "#" === e.value.bondType) ||
                                    ("=" === t.value.bondType &&
                                        e &&
                                        0 === e.value.rings.length &&
                                        "=" === e.value.bondType &&
                                        "-" !== t.value.branchBond)
                                )
                                    ((t.value.drawExplicit = !1),
                                        e && (this.graph.getEdge(t.id, e.id).center = !0),
                                        (this.graph.getEdge(t.id, n.id).center = !0),
                                        ("#" === t.value.bondType || (e && "#" === e.value.bondType)) && (n.angle = 0),
                                        (n.drawExplicit = !0),
                                        this.createNextBond(n, t, u + n.angle));
                                else if (e && e.value.rings.length > 0) {
                                    let e = r.toRad(60),
                                        i = -e,
                                        o = new a(this.opts.bondLength, 0),
                                        s = new a(this.opts.bondLength, 0);
                                    (o.rotate(e).add(t.position), s.rotate(i).add(t.position));
                                    let l = this.getCurrentCenterOfMass(),
                                        c = o.distanceSq(l),
                                        h = s.distanceSq(l);
                                    ((n.angle = c < h ? i : e), this.createNextBond(n, t, u + n.angle));
                                } else {
                                    let r = t.angle;
                                    if (
                                        (e && e.neighbours.length > 3
                                            ? (r = r > 0 ? Math.min(1.0472, r) : r < 0 ? Math.max(-1.0472, r) : 1.0472)
                                            : r || ((r = this.getLastVertexWithAngle(t.id).angle), r || (r = 1.0472)),
                                        e && !l)
                                    ) {
                                        let e = this.graph.getEdge(t.id, n.id).bondType;
                                        "/" === e
                                            ? ("/" === this.doubleBondConfig ||
                                                  ("\\" === this.doubleBondConfig && (r = -r)),
                                              (this.doubleBondConfig = null))
                                            : "\\" === e &&
                                              ("/" === this.doubleBondConfig ? (r = -r) : this.doubleBondConfig,
                                              (this.doubleBondConfig = null));
                                    }
                                    ((n.angle = o ? r : -r), this.createNextBond(n, t, u + n.angle));
                                }
                            } else if (2 === s.length) {
                                let n = t.angle;
                                n || (n = 1.0472);
                                let r = this.graph.getTreeDepth(s[0], t.id),
                                    i = this.graph.getTreeDepth(s[1], t.id),
                                    a = this.graph.vertices[s[0]],
                                    o = this.graph.vertices[s[1]];
                                ((a.value.subtreeDepth = r), (o.value.subtreeDepth = i));
                                let l = this.graph.getTreeDepth(e ? e.id : null, t.id);
                                e && (e.value.subtreeDepth = l);
                                let c = 0,
                                    h = 1;
                                "C" === o.value.element && "C" !== a.value.element && i > 1 && r < 5
                                    ? ((c = 1), (h = 0))
                                    : "C" !== o.value.element && "C" === a.value.element && r > 1 && i < 5
                                      ? ((c = 0), (h = 1))
                                      : i > r && ((c = 1), (h = 0));
                                let d = this.graph.vertices[s[c]],
                                    g = this.graph.vertices[s[h]],
                                    f = (this.graph.getEdge(t.id, d.id), this.graph.getEdge(t.id, g.id), !1);
                                (l < r && l < i && (f = !0),
                                    (g.angle = n),
                                    (d.angle = -n),
                                    "\\" === this.doubleBondConfig
                                        ? "\\" === g.value.branchBond && ((g.angle = -n), (d.angle = n))
                                        : "/" === this.doubleBondConfig &&
                                          "/" === g.value.branchBond &&
                                          ((g.angle = -n), (d.angle = n)),
                                    this.createNextBond(g, t, u + g.angle, f),
                                    this.createNextBond(d, t, u + d.angle, f));
                            } else if (3 === s.length) {
                                let n = this.graph.getTreeDepth(s[0], t.id),
                                    i = this.graph.getTreeDepth(s[1], t.id),
                                    a = this.graph.getTreeDepth(s[2], t.id),
                                    o = this.graph.vertices[s[0]],
                                    l = this.graph.vertices[s[1]],
                                    c = this.graph.vertices[s[2]];
                                ((o.value.subtreeDepth = n),
                                    (l.value.subtreeDepth = i),
                                    (c.value.subtreeDepth = a),
                                    i > n && i > a
                                        ? ((o = this.graph.vertices[s[1]]),
                                          (l = this.graph.vertices[s[0]]),
                                          (c = this.graph.vertices[s[2]]))
                                        : a > n &&
                                          a > i &&
                                          ((o = this.graph.vertices[s[2]]),
                                          (l = this.graph.vertices[s[0]]),
                                          (c = this.graph.vertices[s[1]])),
                                    e &&
                                    e.value.rings.length < 1 &&
                                    o.value.rings.length < 1 &&
                                    l.value.rings.length < 1 &&
                                    c.value.rings.length < 1 &&
                                    1 === this.graph.getTreeDepth(l.id, t.id) &&
                                    1 === this.graph.getTreeDepth(c.id, t.id) &&
                                    this.graph.getTreeDepth(o.id, t.id) > 1
                                        ? ((o.angle = -t.angle),
                                          t.angle >= 0
                                              ? ((l.angle = r.toRad(30)), (c.angle = r.toRad(90)))
                                              : ((l.angle = -r.toRad(30)), (c.angle = -r.toRad(90))),
                                          this.createNextBond(o, t, u + o.angle),
                                          this.createNextBond(l, t, u + l.angle),
                                          this.createNextBond(c, t, u + c.angle))
                                        : ((o.angle = 0),
                                          (l.angle = r.toRad(90)),
                                          (c.angle = -r.toRad(90)),
                                          this.createNextBond(o, t, u + o.angle),
                                          this.createNextBond(l, t, u + l.angle),
                                          this.createNextBond(c, t, u + c.angle)));
                            } else if (4 === s.length) {
                                let e = this.graph.getTreeDepth(s[0], t.id),
                                    n = this.graph.getTreeDepth(s[1], t.id),
                                    i = this.graph.getTreeDepth(s[2], t.id),
                                    a = this.graph.getTreeDepth(s[3], t.id),
                                    o = this.graph.vertices[s[0]],
                                    l = this.graph.vertices[s[1]],
                                    c = this.graph.vertices[s[2]],
                                    h = this.graph.vertices[s[3]];
                                ((o.value.subtreeDepth = e),
                                    (l.value.subtreeDepth = n),
                                    (c.value.subtreeDepth = i),
                                    (h.value.subtreeDepth = a),
                                    n > e && n > i && n > a
                                        ? ((o = this.graph.vertices[s[1]]),
                                          (l = this.graph.vertices[s[0]]),
                                          (c = this.graph.vertices[s[2]]),
                                          (h = this.graph.vertices[s[3]]))
                                        : i > e && i > n && i > a
                                          ? ((o = this.graph.vertices[s[2]]),
                                            (l = this.graph.vertices[s[0]]),
                                            (c = this.graph.vertices[s[1]]),
                                            (h = this.graph.vertices[s[3]]))
                                          : a > e &&
                                            a > n &&
                                            a > i &&
                                            ((o = this.graph.vertices[s[3]]),
                                            (l = this.graph.vertices[s[0]]),
                                            (c = this.graph.vertices[s[1]]),
                                            (h = this.graph.vertices[s[2]])),
                                    (o.angle = -r.toRad(36)),
                                    (l.angle = r.toRad(36)),
                                    (c.angle = -r.toRad(108)),
                                    (h.angle = r.toRad(108)),
                                    this.createNextBond(o, t, u + o.angle),
                                    this.createNextBond(l, t, u + l.angle),
                                    this.createNextBond(c, t, u + c.angle),
                                    this.createNextBond(h, t, u + h.angle));
                            }
                        }
                    }
                    getCommonRingbondNeighbour(t) {
                        let e = t.neighbours;
                        for (var n = 0; n < e.length; n++) {
                            let r = this.graph.vertices[e[n]];
                            if (i.containsAll(r.value.rings, t.value.rings)) return r;
                        }
                        return null;
                    }
                    isPointInRing(t) {
                        for (var e = 0; e < this.rings.length; e++) {
                            let n = this.rings[e];
                            if (!n.positioned) continue;
                            let i = r.polyCircumradius(this.opts.bondLength, n.getSize()),
                                a = i * i;
                            if (t.distanceSq(n.center) < a) return !0;
                        }
                        return !1;
                    }
                    isEdgeInRing(t) {
                        let e = this.graph.vertices[t.sourceId],
                            n = this.graph.vertices[t.targetId];
                        return this.areVerticesInSameRing(e, n);
                    }
                    isEdgeRotatable(t) {
                        let e = this.graph.vertices[t.sourceId],
                            n = this.graph.vertices[t.targetId];
                        return !(
                            "-" !== t.bondType ||
                            e.isTerminal() ||
                            n.isTerminal() ||
                            (e.value.rings.length > 0 && n.value.rings.length > 0 && this.areVerticesInSameRing(e, n))
                        );
                    }
                    isRingAromatic(t) {
                        for (var e = 0; e < t.members.length; e++)
                            if (!this.graph.vertices[t.members[e]].value.isPartOfAromaticRing) return !1;
                        return !0;
                    }
                    getEdgeNormals(t) {
                        let e = this.graph.vertices[t.sourceId].position,
                            n = this.graph.vertices[t.targetId].position;
                        return a.units(e, n);
                    }
                    getNonRingNeighbours(t) {
                        let e = Array(),
                            n = this.graph.vertices[t],
                            r = n.neighbours;
                        for (var a = 0; a < r.length; a++) {
                            let t = this.graph.vertices[r[a]];
                            0 === i.intersection(n.value.rings, t.value.rings).length &&
                                0 == t.value.isBridge &&
                                e.push(t);
                        }
                        return e;
                    }
                    annotateStereochemistry() {
                        for (var t = 0; t < this.graph.vertices.length; t++) {
                            let a = this.graph.vertices[t];
                            if (!a.value.isStereoCenter) continue;
                            let o = a.getNeighbours(),
                                s = o.length,
                                l = Array(s);
                            for (var e = 0; e < s; e++) {
                                let t = new Uint8Array(this.graph.vertices.length),
                                    r = Array(Array());
                                ((t[a.id] = 1), this.visitStereochemistry(o[e], a.id, t, r, 10, 0));
                                for (var n = 0; n < r.length; n++)
                                    r[n].sort(function (t, e) {
                                        return e - t;
                                    });
                                l[e] = [e, r];
                            }
                            let c = 0,
                                u = 0;
                            for (e = 0; e < l.length; e++)
                                for (l[e][1].length > c && (c = l[e][1].length), n = 0; n < l[e][1].length; n++)
                                    l[e][1][n].length > u && (u = l[e][1][n].length);
                            for (e = 0; e < l.length; e++) {
                                let t = c - l[e][1].length;
                                for (n = 0; n < t; n++) l[e][1].push([]);
                                for (l[e][1].push([o[e]]), n = 0; n < l[e][1].length; n++) {
                                    let t = u - l[e][1][n].length;
                                    for (var i = 0; i < t; i++) l[e][1][n].push(0);
                                }
                            }
                            l.sort(function (t, e) {
                                for (var n = 0; n < t[1].length; n++)
                                    for (var r = 0; r < t[1][n].length; r++) {
                                        if (t[1][n][r] > e[1][n][r]) return -1;
                                        if (t[1][n][r] < e[1][n][r]) return 1;
                                    }
                                return 0;
                            });
                            let h = new Uint8Array(s);
                            for (e = 0; e < s; e++) ((h[e] = l[e][0]), (a.value.priority = e));
                            let d = this.graph.vertices[o[h[0]]].position,
                                g = this.graph.vertices[o[h[1]]].position,
                                f = this.graph.vertices[o[h[2]]].position,
                                p = d.relativeClockwise(g, a.position),
                                m = (d.relativeClockwise(f, a.position), -1 === p),
                                v = "@" === a.value.bracket.chirality ? -1 : 1,
                                y = r.parityOfPermutation(h) * v == 1 ? "R" : "S",
                                b = "down",
                                x = "up";
                            (((m && "R" !== y) || (!m && "S" !== y)) &&
                                ((a.value.hydrogenDirection = "up"), (b = "up"), (x = "down")),
                                a.value.hasHydrogen && (this.graph.getEdge(a.id, o[h[h.length - 1]]).wedge = b));
                            let w = new Array(o.length - 1),
                                $ = a.value.rings.length > 1 && a.value.hasHydrogen,
                                A = a.value.hasHydrogen ? 1 : 0;
                            for (e = 0; e < h.length - A; e++) {
                                w[e] = new Uint32Array(2);
                                let t = this.graph.vertices[o[h[e]]];
                                ((w[e][0] += t.value.isStereoCenter ? 0 : 1e5),
                                    (w[e][0] += this.areVerticesInSameRing(t, a) ? 0 : 1e4),
                                    (w[e][0] += t.value.isHeteroAtom() ? 1e3 : 0),
                                    (w[e][0] -= 0 === t.value.subtreeDepth ? 1e3 : 0),
                                    (w[e][0] += 1e3 - t.value.subtreeDepth),
                                    (w[e][1] = o[h[e]]));
                            }
                            if (
                                (w.sort(function (t, e) {
                                    return t[0] > e[0] ? -1 : t[0] < e[0] ? 1 : 0;
                                }),
                                !$)
                            ) {
                                let t = w[0][1];
                                if (a.value.hasHydrogen) this.graph.getEdge(a.id, t).wedge = x;
                                else {
                                    let n = x;
                                    for (e = h.length - 1; e >= 0 && ((n = n === b ? x : b), o[h[e]] !== t); e--);
                                    this.graph.getEdge(a.id, t).wedge = n;
                                }
                            }
                            a.value.chirality = y;
                        }
                    }
                    visitStereochemistry(t, e, n, r, i, a, o = 0) {
                        n[t] = 1;
                        let s = this.graph.vertices[t],
                            l = s.value.getAtomicNumber();
                        r.length <= a && r.push(Array());
                        for (var c = 0; c < this.graph.getEdge(t, e).weight; c++) r[a].push(1e3 * o + l);
                        let u = this.graph.vertices[t].neighbours;
                        for (c = 0; c < u.length; c++)
                            1 !== n[u[c]] && a < i - 1 && this.visitStereochemistry(u[c], t, n.slice(), r, i, a + 1, l);
                        if (a < i - 1) {
                            let e = 0;
                            for (c = 0; c < u.length; c++) e += this.graph.getEdge(t, u[c]).weight;
                            for (c = 0; c < s.value.getMaxBonds() - e; c++)
                                (r.length <= a + 1 && r.push(Array()), r[a + 1].push(1e3 * l + 1));
                        }
                    }
                    initPseudoElements() {
                        for (var t = 0; t < this.graph.vertices.length; t++) {
                            const n = this.graph.vertices[t],
                                r = n.neighbours;
                            let i = Array(r.length);
                            for (var e = 0; e < r.length; e++) i[e] = this.graph.vertices[r[e]];
                            if (n.getNeighbourCount() < 3 || n.value.rings.length > 0) continue;
                            if ("P" === n.value.element) continue;
                            if (
                                "C" === n.value.element &&
                                3 === i.length &&
                                "N" === i[0].value.element &&
                                "N" === i[1].value.element &&
                                "N" === i[2].value.element
                            )
                                continue;
                            let a = 0,
                                o = 0;
                            for (e = 0; e < i.length; e++) {
                                let t = i[e],
                                    n = t.value.element,
                                    r = t.getNeighbourCount();
                                ("C" !== n && "H" !== n && 1 === r && a++, r > 1 && o++);
                            }
                            if (o > 1 || a < 2) continue;
                            let s = null;
                            for (e = 0; e < i.length; e++) {
                                let t = i[e];
                                t.getNeighbourCount() > 1 && (s = t);
                            }
                            for (e = 0; e < i.length; e++) {
                                let t = i[e];
                                if (t.getNeighbourCount() > 1) continue;
                                t.value.isDrawn = !1;
                                let r = l.maxBonds[t.value.element] - t.value.bondCount,
                                    a = "";
                                (t.value.bracket && ((r = t.value.bracket.hcount), (a = t.value.bracket.charge || 0)),
                                    n.value.attachPseudoElement(t.value.element, s ? s.value.element : null, r, a));
                            }
                        }
                        for (t = 0; t < this.graph.vertices.length; t++) {
                            const n = this.graph.vertices[t],
                                r = n.value,
                                i = r.element;
                            if ("C" === i || "H" === i || !r.isDrawn) continue;
                            const a = n.neighbours;
                            let o = Array(a.length);
                            for (e = 0; e < a.length; e++) o[e] = this.graph.vertices[a[e]];
                            for (e = 0; e < o.length; e++) {
                                let t = o[e].value;
                                if (!t.hasAttachedPseudoElements || 2 !== t.getAttachedPseudoElementsCount()) continue;
                                const r = t.getAttachedPseudoElements();
                                r.hasOwnProperty("0O") &&
                                    r.hasOwnProperty("3C") &&
                                    ((t.isDrawn = !1), n.value.attachPseudoElement("Ac", "", 0));
                            }
                        }
                    }
                };
            },
            826: (t) => {
                class e {
                    constructor(t, e, n = 1) {
                        ((this.id = null),
                            (this.sourceId = t),
                            (this.targetId = e),
                            (this.weight = n),
                            (this.bondType = "-"),
                            (this.isPartOfAromaticRing = !1),
                            (this.center = !1),
                            (this.wedge = ""));
                    }
                    setBondType(t) {
                        ((this.bondType = t), (this.weight = e.bonds[t]));
                    }
                    static get bonds() {
                        return { "-": 1, "/": 1, "\\": 1, "=": 2, "#": 3, $: 4 };
                    }
                }
                t.exports = e;
            },
            707: (t, e, n) => {
                const r = n(474),
                    i = (n(614), n(152)),
                    a = n(826),
                    o = (n(421), n(427));
                class s {
                    constructor(t, e = !1) {
                        ((this.vertices = Array()),
                            (this.edges = Array()),
                            (this.vertexIdsToEdgeId = {}),
                            (this.isomeric = e),
                            (this._time = 0),
                            this._init(t));
                    }
                    _init(t, e = 0, n = null, r = !1) {
                        let s = new o(t.atom.element ? t.atom.element : t.atom, t.bond);
                        ((s.branchBond = t.branchBond),
                            (s.ringbonds = t.ringbonds),
                            (s.bracket = t.atom.element ? t.atom : null));
                        let l = new i(s),
                            c = this.vertices[n];
                        if ((this.addVertex(l), null !== n)) {
                            (l.setParentVertexId(n),
                                l.value.addNeighbouringElement(c.value.element),
                                c.addChild(l.id),
                                c.value.addNeighbouringElement(s.element),
                                c.spanningTreeChildren.push(l.id));
                            let t = new a(n, l.id, 1),
                                e = null;
                            (r
                                ? (t.setBondType(l.value.branchBond || "-"),
                                  (e = l.id),
                                  t.setBondType(l.value.branchBond || "-"),
                                  (e = l.id))
                                : (t.setBondType(c.value.bondType || "-"), (e = c.id)),
                                this.addEdge(t));
                        }
                        let u = t.ringbondCount + 1;
                        s.bracket && (u += s.bracket.hcount);
                        let h = 0;
                        if (s.bracket && s.bracket.chirality) {
                            ((s.isStereoCenter = !0), (h = s.bracket.hcount));
                            for (var d = 0; d < h; d++)
                                this._init(
                                    {
                                        atom: "H",
                                        isBracket: "false",
                                        branches: Array(),
                                        branchCount: 0,
                                        ringbonds: Array(),
                                        ringbondCount: !1,
                                        next: null,
                                        hasNext: !1,
                                        bond: "-",
                                    },
                                    d,
                                    l.id,
                                    !0,
                                );
                        }
                        for (d = 0; d < t.branchCount; d++) this._init(t.branches[d], d + u, l.id, !0);
                        t.hasNext && this._init(t.next, t.branchCount + u, l.id);
                    }
                    clear() {
                        ((this.vertices = Array()), (this.edges = Array()), (this.vertexIdsToEdgeId = {}));
                    }
                    addVertex(t) {
                        return ((t.id = this.vertices.length), this.vertices.push(t), t.id);
                    }
                    addEdge(t) {
                        let e = this.vertices[t.sourceId],
                            n = this.vertices[t.targetId];
                        return (
                            (t.id = this.edges.length),
                            this.edges.push(t),
                            (this.vertexIdsToEdgeId[t.sourceId + "_" + t.targetId] = t.id),
                            (this.vertexIdsToEdgeId[t.targetId + "_" + t.sourceId] = t.id),
                            (t.isPartOfAromaticRing = e.value.isPartOfAromaticRing && n.value.isPartOfAromaticRing),
                            (e.value.bondCount += t.weight),
                            (n.value.bondCount += t.weight),
                            e.edges.push(t.id),
                            n.edges.push(t.id),
                            t.id
                        );
                    }
                    getEdge(t, e) {
                        let n = this.vertexIdsToEdgeId[t + "_" + e];
                        return void 0 === n ? null : this.edges[n];
                    }
                    getEdges(t) {
                        let e = Array(),
                            n = this.vertices[t];
                        for (var r = 0; r < n.neighbours.length; r++)
                            e.push(this.vertexIdsToEdgeId[t + "_" + n.neighbours[r]]);
                        return e;
                    }
                    hasEdge(t, e) {
                        return void 0 !== this.vertexIdsToEdgeId[t + "_" + e];
                    }
                    getVertexList() {
                        let t = [this.vertices.length];
                        for (var e = 0; e < this.vertices.length; e++) t[e] = this.vertices[e].id;
                        return t;
                    }
                    getEdgeList() {
                        let t = Array(this.edges.length);
                        for (var e = 0; e < this.edges.length; e++)
                            t[e] = [this.edges[e].sourceId, this.edges[e].targetId];
                        return t;
                    }
                    getAdjacencyMatrix() {
                        let t = this.vertices.length,
                            e = Array(t);
                        for (var n = 0; n < t; n++) ((e[n] = new Array(t)), e[n].fill(0));
                        for (n = 0; n < this.edges.length; n++) {
                            let t = this.edges[n];
                            ((e[t.sourceId][t.targetId] = 1), (e[t.targetId][t.sourceId] = 1));
                        }
                        return e;
                    }
                    getComponentsAdjacencyMatrix() {
                        let t = this.vertices.length,
                            e = Array(t),
                            n = this.getBridges();
                        for (var r = 0; r < t; r++) ((e[r] = new Array(t)), e[r].fill(0));
                        for (r = 0; r < this.edges.length; r++) {
                            let t = this.edges[r];
                            ((e[t.sourceId][t.targetId] = 1), (e[t.targetId][t.sourceId] = 1));
                        }
                        for (r = 0; r < n.length; r++) ((e[n[r][0]][n[r][1]] = 0), (e[n[r][1]][n[r][0]] = 0));
                        return e;
                    }
                    getSubgraphAdjacencyMatrix(t) {
                        let e = t.length,
                            n = Array(e);
                        for (var r = 0; r < e; r++) {
                            ((n[r] = new Array(e)), n[r].fill(0));
                            for (var i = 0; i < e; i++) r !== i && this.hasEdge(t[r], t[i]) && (n[r][i] = 1);
                        }
                        return n;
                    }
                    getDistanceMatrix() {
                        let t = this.vertices.length,
                            e = this.getAdjacencyMatrix(),
                            n = Array(t);
                        for (var r = 0; r < t; r++) ((n[r] = Array(t)), n[r].fill(1 / 0));
                        for (r = 0; r < t; r++) for (var i = 0; i < t; i++) 1 === e[r][i] && (n[r][i] = 1);
                        for (var a = 0; a < t; a++)
                            for (r = 0; r < t; r++)
                                for (i = 0; i < t; i++) n[r][i] > n[r][a] + n[a][i] && (n[r][i] = n[r][a] + n[a][i]);
                        return n;
                    }
                    getSubgraphDistanceMatrix(t) {
                        let e = t.length,
                            n = this.getSubgraphAdjacencyMatrix(t),
                            r = Array(e);
                        for (var i = 0; i < e; i++) ((r[i] = Array(e)), r[i].fill(1 / 0));
                        for (i = 0; i < e; i++) for (var a = 0; a < e; a++) 1 === n[i][a] && (r[i][a] = 1);
                        for (var o = 0; o < e; o++)
                            for (i = 0; i < e; i++)
                                for (a = 0; a < e; a++) r[i][a] > r[i][o] + r[o][a] && (r[i][a] = r[i][o] + r[o][a]);
                        return r;
                    }
                    getAdjacencyList() {
                        let t = this.vertices.length,
                            e = Array(t);
                        for (var n = 0; n < t; n++) {
                            e[n] = [];
                            for (var r = 0; r < t; r++)
                                n !== r && this.hasEdge(this.vertices[n].id, this.vertices[r].id) && e[n].push(r);
                        }
                        return e;
                    }
                    getSubgraphAdjacencyList(t) {
                        let e = t.length,
                            n = Array(e);
                        for (var r = 0; r < e; r++) {
                            n[r] = Array();
                            for (var i = 0; i < e; i++) r !== i && this.hasEdge(t[r], t[i]) && n[r].push(i);
                        }
                        return n;
                    }
                    getBridges() {
                        let t = this.vertices.length,
                            e = new Array(t),
                            n = new Array(t),
                            r = new Array(t),
                            i = new Array(t),
                            a = this.getAdjacencyList(),
                            o = Array();
                        (e.fill(!1), i.fill(null), (this._time = 0));
                        for (var s = 0; s < t; s++) e[s] || this._bridgeDfs(s, e, n, r, i, a, o);
                        return o;
                    }
                    traverseBF(t, e) {
                        let n = this.vertices.length,
                            r = new Array(n);
                        r.fill(!1);
                        for (var i = [t]; i.length > 0; ) {
                            let t = i.shift(),
                                n = this.vertices[t];
                            e(n);
                            for (var a = 0; a < n.neighbours.length; a++) {
                                let t = n.neighbours[a];
                                r[t] || ((r[t] = !0), i.push(t));
                            }
                        }
                    }
                    getTreeDepth(t, e) {
                        if (null === t || null === e) return 0;
                        let n = this.vertices[t].getSpanningTreeNeighbours(e),
                            r = 0;
                        for (var i = 0; i < n.length; i++) {
                            let e = n[i],
                                a = this.getTreeDepth(e, t);
                            a > r && (r = a);
                        }
                        return r + 1;
                    }
                    traverseTree(t, e, n, r = Number.MAX_SAFE_INTEGER, i = !1, a = 1, o = null) {
                        if ((null === o && (o = new Uint8Array(this.vertices.length)), a > r + 1 || 1 === o[t])) return;
                        o[t] = 1;
                        let s = this.vertices[t],
                            l = s.getNeighbours(e);
                        (!i || a > 1) && n(s);
                        for (var c = 0; c < l.length; c++) this.traverseTree(l[c], t, n, r, i, a + 1, o);
                    }
                    kkLayout(t, e, n, i, a) {
                        let o = a;
                        for (var s = t.length; s--; ) var l = this.vertices[t[s]].neighbours.length;
                        let c = this.getSubgraphDistanceMatrix(t),
                            u = t.length,
                            h = r.polyCircumradius(500, u),
                            d = r.centralAngle(u),
                            g = 0,
                            f = new Float32Array(u),
                            p = new Float32Array(u),
                            m = Array(u);
                        for (s = u; s--; ) {
                            let n = this.vertices[t[s]];
                            (n.positioned
                                ? ((f[s] = n.position.x), (p[s] = n.position.y))
                                : ((f[s] = e.x + Math.cos(g) * h), (p[s] = e.y + Math.sin(g) * h)),
                                (m[s] = n.positioned),
                                (g += d));
                        }
                        let v = Array(u);
                        for (s = u; s--; ) for (v[s] = new Array(u), l = u; l--; ) v[s][l] = a * c[s][l];
                        let y = Array(u);
                        for (s = u; s--; ) for (y[s] = Array(u), l = u; l--; ) y[s][l] = o * Math.pow(c[s][l], -2);
                        let b,
                            x,
                            w,
                            $,
                            A,
                            _,
                            C,
                            S = Array(u),
                            M = new Float32Array(u),
                            T = new Float32Array(u);
                        for (s = u; s--; ) S[s] = Array(u);
                        for (s = u; s--; ) {
                            ((b = f[s]), (x = p[s]), (w = 0), ($ = 0));
                            let t = u;
                            for (; t--; )
                                s !== t &&
                                    ((A = f[t]),
                                    (_ = p[t]),
                                    (C = 1 / Math.sqrt((b - A) * (b - A) + (x - _) * (x - _))),
                                    (S[s][t] = [
                                        y[s][t] * (b - A - v[s][t] * (b - A) * C),
                                        y[s][t] * (x - _ - v[s][t] * (x - _) * C),
                                    ]),
                                    (S[t][s] = S[s][t]),
                                    (w += S[s][t][0]),
                                    ($ += S[s][t][1]));
                            ((M[s] = w), (T[s] = $));
                        }
                        let k = function (t) {
                                return [M[t] * M[t] + T[t] * T[t], M[t], T[t]];
                            },
                            R = function () {
                                let t = 0,
                                    e = 0,
                                    n = 0,
                                    r = 0;
                                for (s = u; s--; ) {
                                    let [i, a, o] = k(s);
                                    i > t && !1 === m[s] && ((t = i), (e = s), (n = a), (r = o));
                                }
                                return [e, t, n, r];
                            },
                            N = function (t, e, n) {
                                let r = 0,
                                    i = 0,
                                    a = 0,
                                    o = f[t],
                                    l = p[t],
                                    c = v[t],
                                    h = y[t];
                                for (s = u; s--; ) {
                                    if (s === t) continue;
                                    let e = f[s],
                                        n = p[s],
                                        u = c[s],
                                        d = h[s],
                                        g = (o - e) * (o - e),
                                        m = 1 / Math.pow(g + (l - n) * (l - n), 1.5);
                                    ((r += d * (1 - u * (l - n) * (l - n) * m)),
                                        (i += d * (1 - u * g * m)),
                                        (a += d * (u * (o - e) * (l - n) * m)));
                                }
                                (0 === r && (r = 0.1), 0 === i && (i = 0.1), 0 === a && (a = 0.1));
                                let d = e / r + n / a;
                                d /= a / r - i / a;
                                let g = -(a * d + e) / r;
                                ((f[t] += g), (p[t] += d));
                                let m,
                                    b,
                                    x,
                                    w,
                                    $,
                                    A = S[t];
                                for (e = 0, n = 0, o = f[t], l = p[t], s = u; s--; )
                                    t !== s &&
                                        ((m = f[s]),
                                        (b = p[s]),
                                        (x = A[s][0]),
                                        (w = A[s][1]),
                                        ($ = 1 / Math.sqrt((o - m) * (o - m) + (l - b) * (l - b))),
                                        (g = h[s] * (o - m - c[s] * (o - m) * $)),
                                        (d = h[s] * (l - b - c[s] * (l - b) * $)),
                                        (A[s] = [g, d]),
                                        (e += g),
                                        (n += d),
                                        (M[s] += g - x),
                                        (T[s] += d - w));
                                ((M[t] = e), (T[t] = n));
                            },
                            L = 1e9,
                            P = 0,
                            D = 0,
                            I = 0,
                            E = 0,
                            B = 0,
                            O = 0;
                        for (; L > 0.1 && 2e3 > B; )
                            for (B++, [P, L, D, I] = R(), E = L, O = 0; E > 0.1 && 50 > O; )
                                (O++, N(P, D, I), ([E, D, I] = k(P)));
                        for (s = u; s--; ) {
                            let e = t[s],
                                n = this.vertices[e];
                            ((n.position.x = f[s]),
                                (n.position.y = p[s]),
                                (n.positioned = !0),
                                (n.forcePositioned = !0));
                        }
                    }
                    _bridgeDfs(t, e, n, r, i, a, o) {
                        ((e[t] = !0), (n[t] = r[t] = ++this._time));
                        for (var s = 0; s < a[t].length; s++) {
                            let l = a[t][s];
                            e[l]
                                ? l !== i[t] && (r[t] = Math.min(r[t], n[l]))
                                : ((i[l] = t),
                                  this._bridgeDfs(l, e, n, r, i, a, o),
                                  (r[t] = Math.min(r[t], r[l])),
                                  r[l] > n[t] && o.push([t, l]));
                        }
                    }
                    static getConnectedComponents(t) {
                        let e = t.length,
                            n = new Array(e),
                            r = new Array();
                        n.fill(!1);
                        for (var i = 0; i < e; i++)
                            if (!n[i]) {
                                let e = Array();
                                ((n[i] = !0), e.push(i), s._ccGetDfs(i, n, t, e), e.length > 1 && r.push(e));
                            }
                        return r;
                    }
                    static getConnectedComponentCount(t) {
                        let e = t.length,
                            n = new Array(e),
                            r = 0;
                        n.fill(!1);
                        for (var i = 0; i < e; i++) n[i] || ((n[i] = !0), r++, s._ccCountDfs(i, n, t));
                        return r;
                    }
                    static _ccCountDfs(t, e, n) {
                        for (var r = 0; r < n[t].length; r++)
                            n[t][r] && !e[r] && t !== r && ((e[r] = !0), s._ccCountDfs(r, e, n));
                    }
                    static _ccGetDfs(t, e, n, r) {
                        for (var i = 0; i < n[t].length; i++)
                            n[t][i] && !e[i] && t !== i && ((e[i] = !0), r.push(i), s._ccGetDfs(i, e, n, r));
                    }
                }
                t.exports = s;
            },
            929: (t, e, n) => {
                const r = n(614);
                class i {
                    constructor(t = new r(0, 0), e = new r(0, 0), n = null, i = null, a = !1, o = !1) {
                        ((this.from = t),
                            (this.to = e),
                            (this.elementFrom = n),
                            (this.elementTo = i),
                            (this.chiralFrom = a),
                            (this.chiralTo = o));
                    }
                    clone() {
                        return new i(this.from.clone(), this.to.clone(), this.elementFrom, this.elementTo);
                    }
                    getLength() {
                        return Math.sqrt(Math.pow(this.to.x - this.from.x, 2) + Math.pow(this.to.y - this.from.y, 2));
                    }
                    getAngle() {
                        return r.subtract(this.getRightVector(), this.getLeftVector()).angle();
                    }
                    getRightVector() {
                        return this.from.x < this.to.x ? this.to : this.from;
                    }
                    getLeftVector() {
                        return this.from.x < this.to.x ? this.from : this.to;
                    }
                    getRightElement() {
                        return this.from.x < this.to.x ? this.elementTo : this.elementFrom;
                    }
                    getLeftElement() {
                        return this.from.x < this.to.x ? this.elementFrom : this.elementTo;
                    }
                    getRightChiral() {
                        return this.from.x < this.to.x ? this.chiralTo : this.chiralFrom;
                    }
                    getLeftChiral() {
                        return this.from.x < this.to.x ? this.chiralFrom : this.chiralTo;
                    }
                    setRightVector(t, e) {
                        return (
                            this.from.x < this.to.x
                                ? ((this.to.x = t), (this.to.y = e))
                                : ((this.from.x = t), (this.from.y = e)),
                            this
                        );
                    }
                    setLeftVector(t, e) {
                        return (
                            this.from.x < this.to.x
                                ? ((this.from.x = t), (this.from.y = e))
                                : ((this.to.x = t), (this.to.y = e)),
                            this
                        );
                    }
                    rotateToXAxis() {
                        let t = this.getLeftVector();
                        return (this.setRightVector(t.x + this.getLength(), t.y), this);
                    }
                    rotate(t) {
                        let e = this.getLeftVector(),
                            n = this.getRightVector(),
                            r = Math.sin(t),
                            i = Math.cos(t),
                            a = i * (n.x - e.x) - r * (n.y - e.y) + e.x,
                            o = r * (n.x - e.x) - i * (n.y - e.y) + e.y;
                        return (this.setRightVector(a, o), this);
                    }
                    shortenFrom(t) {
                        let e = r.subtract(this.to, this.from);
                        return (e.normalize(), e.multiplyScalar(t), this.from.add(e), this);
                    }
                    shortenTo(t) {
                        let e = r.subtract(this.from, this.to);
                        return (e.normalize(), e.multiplyScalar(t), this.to.add(e), this);
                    }
                    shortenRight(t) {
                        return (this.from.x < this.to.x ? this.shortenTo(t) : this.shortenFrom(t), this);
                    }
                    shortenLeft(t) {
                        return (this.from.x < this.to.x ? this.shortenFrom(t) : this.shortenTo(t), this);
                    }
                    shorten(t) {
                        let e = r.subtract(this.from, this.to);
                        return (e.normalize(), e.multiplyScalar(t / 2), this.to.add(e), this.from.subtract(e), this);
                    }
                }
                t.exports = i;
            },
            474: (t) => {
                class e {
                    static round(t, e) {
                        return ((e = e || 1), Number(Math.round(t + "e" + e) + "e-" + e));
                    }
                    static meanAngle(t) {
                        let e = 0,
                            n = 0;
                        for (var r = 0; r < t.length; r++) ((e += Math.sin(t[r])), (n += Math.cos(t[r])));
                        return Math.atan2(e / t.length, n / t.length);
                    }
                    static innerAngle(t) {
                        return e.toRad((180 * (t - 2)) / t);
                    }
                    static polyCircumradius(t, e) {
                        return t / (2 * Math.sin(Math.PI / e));
                    }
                    static apothem(t, e) {
                        return t * Math.cos(Math.PI / e);
                    }
                    static apothemFromSideLength(t, n) {
                        let r = e.polyCircumradius(t, n);
                        return e.apothem(r, n);
                    }
                    static centralAngle(t) {
                        return e.toRad(360 / t);
                    }
                    static toDeg(t) {
                        return t * e.degFactor;
                    }
                    static toRad(t) {
                        return t * e.radFactor;
                    }
                    static parityOfPermutation(t) {
                        let e = new Uint8Array(t.length),
                            n = 0,
                            r = function (n, i = 0) {
                                return 1 === e[n] ? i : (i++, (e[n] = 1), r(t[n], i));
                            };
                        for (var i = 0; i < t.length; i++) 1 !== e[i] && (n += 1 - (r(i) % 2));
                        return n % 2 ? -1 : 1;
                    }
                    static get radFactor() {
                        return Math.PI / 180;
                    }
                    static get degFactor() {
                        return 180 / Math.PI;
                    }
                    static get twoPI() {
                        return 2 * Math.PI;
                    }
                }
                t.exports = e;
            },
            19: (t) => {
                t.exports = (function () {
                    "use strict";
                    function t(e, n, r, i) {
                        ((this.message = e),
                            (this.expected = n),
                            (this.found = r),
                            (this.location = i),
                            (this.name = "SyntaxError"),
                            "function" == typeof Error.captureStackTrace && Error.captureStackTrace(this, t));
                    }
                    return (
                        (function (t, e) {
                            function n() {
                                this.constructor = t;
                            }
                            ((n.prototype = e.prototype), (t.prototype = new n()));
                        })(t, Error),
                        (t.buildMessage = function (t, e) {
                            var n = {
                                literal: function (t) {
                                    return '"' + i(t.text) + '"';
                                },
                                class: function (t) {
                                    var e,
                                        n = "";
                                    for (e = 0; e < t.parts.length; e++)
                                        n +=
                                            t.parts[e] instanceof Array
                                                ? a(t.parts[e][0]) + "-" + a(t.parts[e][1])
                                                : a(t.parts[e]);
                                    return "[" + (t.inverted ? "^" : "") + n + "]";
                                },
                                any: function (t) {
                                    return "any character";
                                },
                                end: function (t) {
                                    return "end of input";
                                },
                                other: function (t) {
                                    return t.description;
                                },
                            };
                            function r(t) {
                                return t.charCodeAt(0).toString(16).toUpperCase();
                            }
                            function i(t) {
                                return t
                                    .replace(/\\/g, "\\\\")
                                    .replace(/"/g, '\\"')
                                    .replace(/\0/g, "\\0")
                                    .replace(/\t/g, "\\t")
                                    .replace(/\n/g, "\\n")
                                    .replace(/\r/g, "\\r")
                                    .replace(/[\x00-\x0F]/g, function (t) {
                                        return "\\x0" + r(t);
                                    })
                                    .replace(/[\x10-\x1F\x7F-\x9F]/g, function (t) {
                                        return "\\x" + r(t);
                                    });
                            }
                            function a(t) {
                                return t
                                    .replace(/\\/g, "\\\\")
                                    .replace(/\]/g, "\\]")
                                    .replace(/\^/g, "\\^")
                                    .replace(/-/g, "\\-")
                                    .replace(/\0/g, "\\0")
                                    .replace(/\t/g, "\\t")
                                    .replace(/\n/g, "\\n")
                                    .replace(/\r/g, "\\r")
                                    .replace(/[\x00-\x0F]/g, function (t) {
                                        return "\\x0" + r(t);
                                    })
                                    .replace(/[\x10-\x1F\x7F-\x9F]/g, function (t) {
                                        return "\\x" + r(t);
                                    });
                            }
                            return (
                                "Expected " +
                                (function (t) {
                                    var e,
                                        r,
                                        i,
                                        a = new Array(t.length);
                                    for (e = 0; e < t.length; e++) a[e] = ((i = t[e]), n[i.type](i));
                                    if ((a.sort(), a.length > 0)) {
                                        for (e = 1, r = 1; e < a.length; e++) a[e - 1] !== a[e] && ((a[r] = a[e]), r++);
                                        a.length = r;
                                    }
                                    switch (a.length) {
                                        case 1:
                                            return a[0];
                                        case 2:
                                            return a[0] + " or " + a[1];
                                        default:
                                            return a.slice(0, -1).join(", ") + ", or " + a[a.length - 1];
                                    }
                                })(t) +
                                " but " +
                                (function (t) {
                                    return t ? '"' + i(t) + '"' : "end of input";
                                })(e) +
                                " found."
                            );
                        }),
                        {
                            SyntaxError: t,
                            parse: function (e, n) {
                                if (((n = void 0 !== n ? n : {}), e.split("(").length - 1 != e.split(")").length - 1))
                                    throw new t(
                                        "The number of opening parentheses does not match the number of closing parentheses.",
                                        null,
                                        null,
                                        0,
                                    );
                                var r,
                                    i = {},
                                    a = { chain: jt },
                                    o = jt,
                                    s = function (t) {
                                        for (var e = [], n = [], r = 0; r < t[1].length; r++) e.push(t[1][r]);
                                        for (r = 0; r < t[2].length; r++) {
                                            var i = t[2][r][0] ? t[2][r][0] : "-";
                                            n.push({ bond: i, id: t[2][r][1] });
                                        }
                                        for (r = 0; r < t[3].length; r++) e.push(t[3][r]);
                                        for (r = 0; r < t[6].length; r++) e.push(t[6][r]);
                                        return {
                                            atom: t[0],
                                            isBracket: !!t[0].element,
                                            branches: e,
                                            branchCount: e.length,
                                            ringbonds: n,
                                            ringbondCount: n.length,
                                            bond: t[4] ? t[4] : "-",
                                            next: t[5],
                                            hasNext: !!t[5],
                                        };
                                    },
                                    l = "(",
                                    c = Ot("(", !1),
                                    u = ")",
                                    h = Ot(")", !1),
                                    d = function (t) {
                                        var e = t[1] ? t[1] : "-";
                                        return ((t[2].branchBond = e), t[2]);
                                    },
                                    g = function (t) {
                                        return t;
                                    },
                                    f = /^[\-=#$:\/\\.]/,
                                    p = Ft(["-", "=", "#", "$", ":", "/", "\\", "."], !1, !1),
                                    m = function (t) {
                                        return t;
                                    },
                                    v = "[",
                                    y = Ot("[", !1),
                                    b = "se",
                                    x = Ot("se", !1),
                                    w = "as",
                                    $ = Ot("as", !1),
                                    A = "]",
                                    _ = Ot("]", !1),
                                    C = function (t) {
                                        return {
                                            isotope: t[1],
                                            element: t[2],
                                            chirality: t[3],
                                            hcount: t[4],
                                            charge: t[5],
                                            class: t[6],
                                        };
                                    },
                                    S = "B",
                                    M = Ot("B", !1),
                                    T = "r",
                                    k = Ot("r", !1),
                                    R = "C",
                                    N = Ot("C", !1),
                                    L = "l",
                                    P = Ot("l", !1),
                                    D = /^[NOPSFI]/,
                                    I = Ft(["N", "O", "P", "S", "F", "I"], !1, !1),
                                    E = function (t) {
                                        return t.length > 1 ? t.join("") : t;
                                    },
                                    B = /^[bcnops]/,
                                    O = Ft(["b", "c", "n", "o", "p", "s"], !1, !1),
                                    F = "*",
                                    z = Ot("*", !1),
                                    q = function (t) {
                                        return t;
                                    },
                                    H = /^[A-Z]/,
                                    j = Ft([["A", "Z"]], !1, !1),
                                    U = /^[a-z]/,
                                    V = Ft([["a", "z"]], !1, !1),
                                    W = function (t) {
                                        return t.join("");
                                    },
                                    Y = "%",
                                    Z = Ot("%", !1),
                                    X = /^[1-9]/,
                                    G = Ft([["1", "9"]], !1, !1),
                                    K = /^[0-9]/,
                                    Q = Ft([["0", "9"]], !1, !1),
                                    J = function (t) {
                                        return 1 == t.length ? Number(t) : Number(t.join("").replace("%", ""));
                                    },
                                    tt = "@",
                                    et = Ot("@", !1),
                                    nt = "TH",
                                    rt = Ot("TH", !1),
                                    it = /^[12]/,
                                    at = Ft(["1", "2"], !1, !1),
                                    ot = "AL",
                                    st = Ot("AL", !1),
                                    lt = "SP",
                                    ct = Ot("SP", !1),
                                    ut = /^[1-3]/,
                                    ht = Ft([["1", "3"]], !1, !1),
                                    dt = "TB",
                                    gt = Ot("TB", !1),
                                    ft = "OH",
                                    pt = Ot("OH", !1),
                                    mt = function (t) {
                                        return t[1] ? ("@" == t[1] ? "@@" : t[1].join("").replace(",", "")) : "@";
                                    },
                                    vt = function (t) {
                                        return t;
                                    },
                                    yt = "+",
                                    bt = Ot("+", !1),
                                    xt = function (t) {
                                        return t[1] ? ("+" != t[1] ? Number(t[1].join("")) : 2) : 1;
                                    },
                                    wt = "-",
                                    $t = Ot("-", !1),
                                    At = function (t) {
                                        return t[1] ? ("-" != t[1] ? -Number(t[1].join("")) : -2) : -1;
                                    },
                                    _t = "H",
                                    Ct = Ot("H", !1),
                                    St = function (t) {
                                        return t[1] ? Number(t[1]) : 1;
                                    },
                                    Mt = ":",
                                    Tt = Ot(":", !1),
                                    kt = /^[0]/,
                                    Rt = Ft(["0"], !1, !1),
                                    Nt = function (t) {
                                        return Number(t[1][0] + t[1][1].join(""));
                                    },
                                    Lt = function (t) {
                                        return Number(t.join(""));
                                    },
                                    Pt = 0,
                                    Dt = [{ line: 1, column: 1 }],
                                    It = 0,
                                    Et = [],
                                    Bt = 0;
                                if ("startRule" in n) {
                                    if (!(n.startRule in a))
                                        throw new Error("Can't start parsing from rule \"" + n.startRule + '".');
                                    o = a[n.startRule];
                                }
                                function Ot(t, e) {
                                    return { type: "literal", text: t, ignoreCase: e };
                                }
                                function Ft(t, e, n) {
                                    return { type: "class", parts: t, inverted: e, ignoreCase: n };
                                }
                                function zt(t) {
                                    var n,
                                        r = Dt[t];
                                    if (r) return r;
                                    for (n = t - 1; !Dt[n]; ) n--;
                                    for (r = { line: (r = Dt[n]).line, column: r.column }; n < t; )
                                        (10 === e.charCodeAt(n) ? (r.line++, (r.column = 1)) : r.column++, n++);
                                    return ((Dt[t] = r), r);
                                }
                                function qt(t, e) {
                                    var n = zt(t),
                                        r = zt(e);
                                    return {
                                        start: { offset: t, line: n.line, column: n.column },
                                        end: { offset: e, line: r.line, column: r.column },
                                    };
                                }
                                function Ht(t) {
                                    Pt < It || (Pt > It && ((It = Pt), (Et = [])), Et.push(t));
                                }
                                function jt() {
                                    var t, n, r, a, o, l, c, u, h;
                                    if (
                                        (Pt,
                                        (t = Pt),
                                        (n = (function () {
                                            var t;
                                            return (
                                                Pt,
                                                (t = (function () {
                                                    var t, n, r;
                                                    return (
                                                        Pt,
                                                        (t = Pt),
                                                        66 === e.charCodeAt(Pt)
                                                            ? ((n = S), Pt++)
                                                            : ((n = i), 0 === Bt && Ht(M)),
                                                        n !== i
                                                            ? (114 === e.charCodeAt(Pt)
                                                                  ? ((r = T), Pt++)
                                                                  : ((r = i), 0 === Bt && Ht(k)),
                                                              r === i && (r = null),
                                                              r !== i ? (t = n = [n, r]) : ((Pt = t), (t = i)))
                                                            : ((Pt = t), (t = i)),
                                                        t === i &&
                                                            ((t = Pt),
                                                            67 === e.charCodeAt(Pt)
                                                                ? ((n = R), Pt++)
                                                                : ((n = i), 0 === Bt && Ht(N)),
                                                            n !== i
                                                                ? (108 === e.charCodeAt(Pt)
                                                                      ? ((r = L), Pt++)
                                                                      : ((r = i), 0 === Bt && Ht(P)),
                                                                  r === i && (r = null),
                                                                  r !== i ? (t = n = [n, r]) : ((Pt = t), (t = i)))
                                                                : ((Pt = t), (t = i)),
                                                            t === i &&
                                                                (D.test(e.charAt(Pt))
                                                                    ? ((t = e.charAt(Pt)), Pt++)
                                                                    : ((t = i), 0 === Bt && Ht(I)))),
                                                        t !== i && (t = E(t)),
                                                        t
                                                    );
                                                })()),
                                                t === i &&
                                                    (t = Wt()) === i &&
                                                    ((t = (function () {
                                                        var t, n, r, a, o, s, l, c, u;
                                                        return (
                                                            Pt,
                                                            (t = Pt),
                                                            91 === e.charCodeAt(Pt)
                                                                ? ((n = v), Pt++)
                                                                : ((n = i), 0 === Bt && Ht(y)),
                                                            n !== i
                                                                ? ((r = (function () {
                                                                      var t, n, r, a;
                                                                      return (
                                                                          Pt,
                                                                          (t = Pt),
                                                                          X.test(e.charAt(Pt))
                                                                              ? ((n = e.charAt(Pt)), Pt++)
                                                                              : ((n = i), 0 === Bt && Ht(G)),
                                                                          n !== i
                                                                              ? (K.test(e.charAt(Pt))
                                                                                    ? ((r = e.charAt(Pt)), Pt++)
                                                                                    : ((r = i), 0 === Bt && Ht(Q)),
                                                                                r === i && (r = null),
                                                                                r !== i
                                                                                    ? (K.test(e.charAt(Pt))
                                                                                          ? ((a = e.charAt(Pt)), Pt++)
                                                                                          : ((a = i),
                                                                                            0 === Bt && Ht(Q)),
                                                                                      a === i && (a = null),
                                                                                      a !== i
                                                                                          ? (t = n = [n, r, a])
                                                                                          : ((Pt = t), (t = i)))
                                                                                    : ((Pt = t), (t = i)))
                                                                              : ((Pt = t), (t = i)),
                                                                          t !== i && (t = Lt(t)),
                                                                          t
                                                                      );
                                                                  })()),
                                                                  r === i && (r = null),
                                                                  r !== i
                                                                      ? (e.substr(Pt, 2) === b
                                                                            ? ((a = b), (Pt += 2))
                                                                            : ((a = i), 0 === Bt && Ht(x)),
                                                                        a === i &&
                                                                            (e.substr(Pt, 2) === w
                                                                                ? ((a = w), (Pt += 2))
                                                                                : ((a = i), 0 === Bt && Ht($)),
                                                                            a === i &&
                                                                                (a = Wt()) === i &&
                                                                                ((a = (function () {
                                                                                    var t, n, r;
                                                                                    return (
                                                                                        Pt,
                                                                                        (t = Pt),
                                                                                        H.test(e.charAt(Pt))
                                                                                            ? ((n = e.charAt(Pt)), Pt++)
                                                                                            : ((n = i),
                                                                                              0 === Bt && Ht(j)),
                                                                                        n !== i
                                                                                            ? (U.test(e.charAt(Pt))
                                                                                                  ? ((r = e.charAt(Pt)),
                                                                                                    Pt++)
                                                                                                  : ((r = i),
                                                                                                    0 === Bt && Ht(V)),
                                                                                              r === i && (r = null),
                                                                                              r !== i
                                                                                                  ? (t = n = [n, r])
                                                                                                  : ((Pt = t), (t = i)))
                                                                                            : ((Pt = t), (t = i)),
                                                                                        t !== i && (t = W(t)),
                                                                                        t
                                                                                    );
                                                                                })()),
                                                                                a === i && (a = Yt()))),
                                                                        a !== i
                                                                            ? ((o = (function () {
                                                                                  var t, n, r, a, o, s;
                                                                                  return (
                                                                                      Pt,
                                                                                      (t = Pt),
                                                                                      64 === e.charCodeAt(Pt)
                                                                                          ? ((n = tt), Pt++)
                                                                                          : ((n = i),
                                                                                            0 === Bt && Ht(et)),
                                                                                      n !== i
                                                                                          ? (64 === e.charCodeAt(Pt)
                                                                                                ? ((r = tt), Pt++)
                                                                                                : ((r = i),
                                                                                                  0 === Bt && Ht(et)),
                                                                                            r === i &&
                                                                                                ((r = Pt),
                                                                                                e.substr(Pt, 2) === nt
                                                                                                    ? ((a = nt),
                                                                                                      (Pt += 2))
                                                                                                    : ((a = i),
                                                                                                      0 === Bt &&
                                                                                                          Ht(rt)),
                                                                                                a !== i
                                                                                                    ? (it.test(
                                                                                                          e.charAt(Pt),
                                                                                                      )
                                                                                                          ? ((o =
                                                                                                                e.charAt(
                                                                                                                    Pt,
                                                                                                                )),
                                                                                                            Pt++)
                                                                                                          : ((o = i),
                                                                                                            0 === Bt &&
                                                                                                                Ht(at)),
                                                                                                      o !== i
                                                                                                          ? (r = a =
                                                                                                                [a, o])
                                                                                                          : ((Pt = r),
                                                                                                            (r = i)))
                                                                                                    : ((Pt = r),
                                                                                                      (r = i)),
                                                                                                r === i &&
                                                                                                    ((r = Pt),
                                                                                                    e.substr(Pt, 2) ===
                                                                                                    ot
                                                                                                        ? ((a = ot),
                                                                                                          (Pt += 2))
                                                                                                        : ((a = i),
                                                                                                          0 === Bt &&
                                                                                                              Ht(st)),
                                                                                                    a !== i
                                                                                                        ? (it.test(
                                                                                                              e.charAt(
                                                                                                                  Pt,
                                                                                                              ),
                                                                                                          )
                                                                                                              ? ((o =
                                                                                                                    e.charAt(
                                                                                                                        Pt,
                                                                                                                    )),
                                                                                                                Pt++)
                                                                                                              : ((o =
                                                                                                                    i),
                                                                                                                0 ===
                                                                                                                    Bt &&
                                                                                                                    Ht(
                                                                                                                        at,
                                                                                                                    )),
                                                                                                          o !== i
                                                                                                              ? (r = a =
                                                                                                                    [
                                                                                                                        a,
                                                                                                                        o,
                                                                                                                    ])
                                                                                                              : ((Pt =
                                                                                                                    r),
                                                                                                                (r =
                                                                                                                    i)))
                                                                                                        : ((Pt = r),
                                                                                                          (r = i)),
                                                                                                    r === i &&
                                                                                                        ((r = Pt),
                                                                                                        e.substr(
                                                                                                            Pt,
                                                                                                            2,
                                                                                                        ) === lt
                                                                                                            ? ((a = lt),
                                                                                                              (Pt += 2))
                                                                                                            : ((a = i),
                                                                                                              0 ===
                                                                                                                  Bt &&
                                                                                                                  Ht(
                                                                                                                      ct,
                                                                                                                  )),
                                                                                                        a !== i
                                                                                                            ? (ut.test(
                                                                                                                  e.charAt(
                                                                                                                      Pt,
                                                                                                                  ),
                                                                                                              )
                                                                                                                  ? ((o =
                                                                                                                        e.charAt(
                                                                                                                            Pt,
                                                                                                                        )),
                                                                                                                    Pt++)
                                                                                                                  : ((o =
                                                                                                                        i),
                                                                                                                    0 ===
                                                                                                                        Bt &&
                                                                                                                        Ht(
                                                                                                                            ht,
                                                                                                                        )),
                                                                                                              o !== i
                                                                                                                  ? (r =
                                                                                                                        a =
                                                                                                                            [
                                                                                                                                a,
                                                                                                                                o,
                                                                                                                            ])
                                                                                                                  : ((Pt =
                                                                                                                        r),
                                                                                                                    (r =
                                                                                                                        i)))
                                                                                                            : ((Pt = r),
                                                                                                              (r = i)),
                                                                                                        r === i &&
                                                                                                            ((r = Pt),
                                                                                                            e.substr(
                                                                                                                Pt,
                                                                                                                2,
                                                                                                            ) === dt
                                                                                                                ? ((a =
                                                                                                                      dt),
                                                                                                                  (Pt += 2))
                                                                                                                : ((a =
                                                                                                                      i),
                                                                                                                  0 ===
                                                                                                                      Bt &&
                                                                                                                      Ht(
                                                                                                                          gt,
                                                                                                                      )),
                                                                                                            a !== i
                                                                                                                ? (X.test(
                                                                                                                      e.charAt(
                                                                                                                          Pt,
                                                                                                                      ),
                                                                                                                  )
                                                                                                                      ? ((o =
                                                                                                                            e.charAt(
                                                                                                                                Pt,
                                                                                                                            )),
                                                                                                                        Pt++)
                                                                                                                      : ((o =
                                                                                                                            i),
                                                                                                                        0 ===
                                                                                                                            Bt &&
                                                                                                                            Ht(
                                                                                                                                G,
                                                                                                                            )),
                                                                                                                  o !==
                                                                                                                  i
                                                                                                                      ? (K.test(
                                                                                                                            e.charAt(
                                                                                                                                Pt,
                                                                                                                            ),
                                                                                                                        )
                                                                                                                            ? ((s =
                                                                                                                                  e.charAt(
                                                                                                                                      Pt,
                                                                                                                                  )),
                                                                                                                              Pt++)
                                                                                                                            : ((s =
                                                                                                                                  i),
                                                                                                                              0 ===
                                                                                                                                  Bt &&
                                                                                                                                  Ht(
                                                                                                                                      Q,
                                                                                                                                  )),
                                                                                                                        s ===
                                                                                                                            i &&
                                                                                                                            (s =
                                                                                                                                null),
                                                                                                                        s !==
                                                                                                                        i
                                                                                                                            ? (r =
                                                                                                                                  a =
                                                                                                                                      [
                                                                                                                                          a,
                                                                                                                                          o,
                                                                                                                                          s,
                                                                                                                                      ])
                                                                                                                            : ((Pt =
                                                                                                                                  r),
                                                                                                                              (r =
                                                                                                                                  i)))
                                                                                                                      : ((Pt =
                                                                                                                            r),
                                                                                                                        (r =
                                                                                                                            i)))
                                                                                                                : ((Pt =
                                                                                                                      r),
                                                                                                                  (r =
                                                                                                                      i)),
                                                                                                            r === i &&
                                                                                                                ((r =
                                                                                                                    Pt),
                                                                                                                e.substr(
                                                                                                                    Pt,
                                                                                                                    2,
                                                                                                                ) === ft
                                                                                                                    ? ((a =
                                                                                                                          ft),
                                                                                                                      (Pt += 2))
                                                                                                                    : ((a =
                                                                                                                          i),
                                                                                                                      0 ===
                                                                                                                          Bt &&
                                                                                                                          Ht(
                                                                                                                              pt,
                                                                                                                          )),
                                                                                                                a !== i
                                                                                                                    ? (X.test(
                                                                                                                          e.charAt(
                                                                                                                              Pt,
                                                                                                                          ),
                                                                                                                      )
                                                                                                                          ? ((o =
                                                                                                                                e.charAt(
                                                                                                                                    Pt,
                                                                                                                                )),
                                                                                                                            Pt++)
                                                                                                                          : ((o =
                                                                                                                                i),
                                                                                                                            0 ===
                                                                                                                                Bt &&
                                                                                                                                Ht(
                                                                                                                                    G,
                                                                                                                                )),
                                                                                                                      o !==
                                                                                                                      i
                                                                                                                          ? (K.test(
                                                                                                                                e.charAt(
                                                                                                                                    Pt,
                                                                                                                                ),
                                                                                                                            )
                                                                                                                                ? ((s =
                                                                                                                                      e.charAt(
                                                                                                                                          Pt,
                                                                                                                                      )),
                                                                                                                                  Pt++)
                                                                                                                                : ((s =
                                                                                                                                      i),
                                                                                                                                  0 ===
                                                                                                                                      Bt &&
                                                                                                                                      Ht(
                                                                                                                                          Q,
                                                                                                                                      )),
                                                                                                                            s ===
                                                                                                                                i &&
                                                                                                                                (s =
                                                                                                                                    null),
                                                                                                                            s !==
                                                                                                                            i
                                                                                                                                ? (r =
                                                                                                                                      a =
                                                                                                                                          [
                                                                                                                                              a,
                                                                                                                                              o,
                                                                                                                                              s,
                                                                                                                                          ])
                                                                                                                                : ((Pt =
                                                                                                                                      r),
                                                                                                                                  (r =
                                                                                                                                      i)))
                                                                                                                          : ((Pt =
                                                                                                                                r),
                                                                                                                            (r =
                                                                                                                                i)))
                                                                                                                    : ((Pt =
                                                                                                                          r),
                                                                                                                      (r =
                                                                                                                          i))))))),
                                                                                            r === i && (r = null),
                                                                                            r !== i
                                                                                                ? (t = n = [n, r])
                                                                                                : ((Pt = t), (t = i)))
                                                                                          : ((Pt = t), (t = i)),
                                                                                      t !== i && (t = mt(t)),
                                                                                      t
                                                                                  );
                                                                              })()),
                                                                              o === i && (o = null),
                                                                              o !== i
                                                                                  ? ((s = (function () {
                                                                                        var t, n, r;
                                                                                        return (
                                                                                            Pt,
                                                                                            (t = Pt),
                                                                                            72 === e.charCodeAt(Pt)
                                                                                                ? ((n = _t), Pt++)
                                                                                                : ((n = i),
                                                                                                  0 === Bt && Ht(Ct)),
                                                                                            n !== i
                                                                                                ? (K.test(e.charAt(Pt))
                                                                                                      ? ((r =
                                                                                                            e.charAt(
                                                                                                                Pt,
                                                                                                            )),
                                                                                                        Pt++)
                                                                                                      : ((r = i),
                                                                                                        0 === Bt &&
                                                                                                            Ht(Q)),
                                                                                                  r === i && (r = null),
                                                                                                  r !== i
                                                                                                      ? (t = n = [n, r])
                                                                                                      : ((Pt = t),
                                                                                                        (t = i)))
                                                                                                : ((Pt = t), (t = i)),
                                                                                            t !== i && (t = St(t)),
                                                                                            t
                                                                                        );
                                                                                    })()),
                                                                                    s === i && (s = null),
                                                                                    s !== i
                                                                                        ? ((l = (function () {
                                                                                              var t;
                                                                                              return (
                                                                                                  Pt,
                                                                                                  (t = (function () {
                                                                                                      var t, n, r, a, o;
                                                                                                      return (
                                                                                                          Pt,
                                                                                                          (t = Pt),
                                                                                                          43 ===
                                                                                                          e.charCodeAt(
                                                                                                              Pt,
                                                                                                          )
                                                                                                              ? ((n =
                                                                                                                    yt),
                                                                                                                Pt++)
                                                                                                              : ((n =
                                                                                                                    i),
                                                                                                                0 ===
                                                                                                                    Bt &&
                                                                                                                    Ht(
                                                                                                                        bt,
                                                                                                                    )),
                                                                                                          n !== i
                                                                                                              ? (43 ===
                                                                                                                e.charCodeAt(
                                                                                                                    Pt,
                                                                                                                )
                                                                                                                    ? ((r =
                                                                                                                          yt),
                                                                                                                      Pt++)
                                                                                                                    : ((r =
                                                                                                                          i),
                                                                                                                      0 ===
                                                                                                                          Bt &&
                                                                                                                          Ht(
                                                                                                                              bt,
                                                                                                                          )),
                                                                                                                r ===
                                                                                                                    i &&
                                                                                                                    ((r =
                                                                                                                        Pt),
                                                                                                                    X.test(
                                                                                                                        e.charAt(
                                                                                                                            Pt,
                                                                                                                        ),
                                                                                                                    )
                                                                                                                        ? ((a =
                                                                                                                              e.charAt(
                                                                                                                                  Pt,
                                                                                                                              )),
                                                                                                                          Pt++)
                                                                                                                        : ((a =
                                                                                                                              i),
                                                                                                                          0 ===
                                                                                                                              Bt &&
                                                                                                                              Ht(
                                                                                                                                  G,
                                                                                                                              )),
                                                                                                                    a !==
                                                                                                                    i
                                                                                                                        ? (K.test(
                                                                                                                              e.charAt(
                                                                                                                                  Pt,
                                                                                                                              ),
                                                                                                                          )
                                                                                                                              ? ((o =
                                                                                                                                    e.charAt(
                                                                                                                                        Pt,
                                                                                                                                    )),
                                                                                                                                Pt++)
                                                                                                                              : ((o =
                                                                                                                                    i),
                                                                                                                                0 ===
                                                                                                                                    Bt &&
                                                                                                                                    Ht(
                                                                                                                                        Q,
                                                                                                                                    )),
                                                                                                                          o ===
                                                                                                                              i &&
                                                                                                                              (o =
                                                                                                                                  null),
                                                                                                                          o !==
                                                                                                                          i
                                                                                                                              ? (r =
                                                                                                                                    a =
                                                                                                                                        [
                                                                                                                                            a,
                                                                                                                                            o,
                                                                                                                                        ])
                                                                                                                              : ((Pt =
                                                                                                                                    r),
                                                                                                                                (r =
                                                                                                                                    i)))
                                                                                                                        : ((Pt =
                                                                                                                              r),
                                                                                                                          (r =
                                                                                                                              i))),
                                                                                                                r ===
                                                                                                                    i &&
                                                                                                                    (r =
                                                                                                                        null),
                                                                                                                r !== i
                                                                                                                    ? (t =
                                                                                                                          n =
                                                                                                                              [
                                                                                                                                  n,
                                                                                                                                  r,
                                                                                                                              ])
                                                                                                                    : ((Pt =
                                                                                                                          t),
                                                                                                                      (t =
                                                                                                                          i)))
                                                                                                              : ((Pt =
                                                                                                                    t),
                                                                                                                (t =
                                                                                                                    i)),
                                                                                                          t !== i &&
                                                                                                              (t =
                                                                                                                  xt(
                                                                                                                      t,
                                                                                                                  )),
                                                                                                          t
                                                                                                      );
                                                                                                  })()),
                                                                                                  t === i &&
                                                                                                      (t =
                                                                                                          (function () {
                                                                                                              var t,
                                                                                                                  n,
                                                                                                                  r,
                                                                                                                  a,
                                                                                                                  o;
                                                                                                              return (
                                                                                                                  Pt,
                                                                                                                  (t =
                                                                                                                      Pt),
                                                                                                                  45 ===
                                                                                                                  e.charCodeAt(
                                                                                                                      Pt,
                                                                                                                  )
                                                                                                                      ? ((n =
                                                                                                                            wt),
                                                                                                                        Pt++)
                                                                                                                      : ((n =
                                                                                                                            i),
                                                                                                                        0 ===
                                                                                                                            Bt &&
                                                                                                                            Ht(
                                                                                                                                $t,
                                                                                                                            )),
                                                                                                                  n !==
                                                                                                                  i
                                                                                                                      ? (45 ===
                                                                                                                        e.charCodeAt(
                                                                                                                            Pt,
                                                                                                                        )
                                                                                                                            ? ((r =
                                                                                                                                  wt),
                                                                                                                              Pt++)
                                                                                                                            : ((r =
                                                                                                                                  i),
                                                                                                                              0 ===
                                                                                                                                  Bt &&
                                                                                                                                  Ht(
                                                                                                                                      $t,
                                                                                                                                  )),
                                                                                                                        r ===
                                                                                                                            i &&
                                                                                                                            ((r =
                                                                                                                                Pt),
                                                                                                                            X.test(
                                                                                                                                e.charAt(
                                                                                                                                    Pt,
                                                                                                                                ),
                                                                                                                            )
                                                                                                                                ? ((a =
                                                                                                                                      e.charAt(
                                                                                                                                          Pt,
                                                                                                                                      )),
                                                                                                                                  Pt++)
                                                                                                                                : ((a =
                                                                                                                                      i),
                                                                                                                                  0 ===
                                                                                                                                      Bt &&
                                                                                                                                      Ht(
                                                                                                                                          G,
                                                                                                                                      )),
                                                                                                                            a !==
                                                                                                                            i
                                                                                                                                ? (K.test(
                                                                                                                                      e.charAt(
                                                                                                                                          Pt,
                                                                                                                                      ),
                                                                                                                                  )
                                                                                                                                      ? ((o =
                                                                                                                                            e.charAt(
                                                                                                                                                Pt,
                                                                                                                                            )),
                                                                                                                                        Pt++)
                                                                                                                                      : ((o =
                                                                                                                                            i),
                                                                                                                                        0 ===
                                                                                                                                            Bt &&
                                                                                                                                            Ht(
                                                                                                                                                Q,
                                                                                                                                            )),
                                                                                                                                  o ===
                                                                                                                                      i &&
                                                                                                                                      (o =
                                                                                                                                          null),
                                                                                                                                  o !==
                                                                                                                                  i
                                                                                                                                      ? (r =
                                                                                                                                            a =
                                                                                                                                                [
                                                                                                                                                    a,
                                                                                                                                                    o,
                                                                                                                                                ])
                                                                                                                                      : ((Pt =
                                                                                                                                            r),
                                                                                                                                        (r =
                                                                                                                                            i)))
                                                                                                                                : ((Pt =
                                                                                                                                      r),
                                                                                                                                  (r =
                                                                                                                                      i))),
                                                                                                                        r ===
                                                                                                                            i &&
                                                                                                                            (r =
                                                                                                                                null),
                                                                                                                        r !==
                                                                                                                        i
                                                                                                                            ? (t =
                                                                                                                                  n =
                                                                                                                                      [
                                                                                                                                          n,
                                                                                                                                          r,
                                                                                                                                      ])
                                                                                                                            : ((Pt =
                                                                                                                                  t),
                                                                                                                              (t =
                                                                                                                                  i)))
                                                                                                                      : ((Pt =
                                                                                                                            t),
                                                                                                                        (t =
                                                                                                                            i)),
                                                                                                                  t !==
                                                                                                                      i &&
                                                                                                                      (t =
                                                                                                                          At(
                                                                                                                              t,
                                                                                                                          )),
                                                                                                                  t
                                                                                                              );
                                                                                                          })()),
                                                                                                  t !== i &&
                                                                                                      (t = vt(t)),
                                                                                                  t
                                                                                              );
                                                                                          })()),
                                                                                          l === i && (l = null),
                                                                                          l !== i
                                                                                              ? ((c = (function () {
                                                                                                    var t,
                                                                                                        n,
                                                                                                        r,
                                                                                                        a,
                                                                                                        o,
                                                                                                        s;
                                                                                                    if (
                                                                                                        (Pt,
                                                                                                        (t = Pt),
                                                                                                        58 ===
                                                                                                        e.charCodeAt(Pt)
                                                                                                            ? ((n = Mt),
                                                                                                              Pt++)
                                                                                                            : ((n = i),
                                                                                                              0 ===
                                                                                                                  Bt &&
                                                                                                                  Ht(
                                                                                                                      Tt,
                                                                                                                  )),
                                                                                                        n !== i)
                                                                                                    ) {
                                                                                                        if (
                                                                                                            ((r = Pt),
                                                                                                            X.test(
                                                                                                                e.charAt(
                                                                                                                    Pt,
                                                                                                                ),
                                                                                                            )
                                                                                                                ? ((a =
                                                                                                                      e.charAt(
                                                                                                                          Pt,
                                                                                                                      )),
                                                                                                                  Pt++)
                                                                                                                : ((a =
                                                                                                                      i),
                                                                                                                  0 ===
                                                                                                                      Bt &&
                                                                                                                      Ht(
                                                                                                                          G,
                                                                                                                      )),
                                                                                                            a !== i)
                                                                                                        ) {
                                                                                                            for (
                                                                                                                o = [],
                                                                                                                    K.test(
                                                                                                                        e.charAt(
                                                                                                                            Pt,
                                                                                                                        ),
                                                                                                                    )
                                                                                                                        ? ((s =
                                                                                                                              e.charAt(
                                                                                                                                  Pt,
                                                                                                                              )),
                                                                                                                          Pt++)
                                                                                                                        : ((s =
                                                                                                                              i),
                                                                                                                          0 ===
                                                                                                                              Bt &&
                                                                                                                              Ht(
                                                                                                                                  Q,
                                                                                                                              ));
                                                                                                                s !== i;
                                                                                                            )
                                                                                                                (o.push(
                                                                                                                    s,
                                                                                                                ),
                                                                                                                    K.test(
                                                                                                                        e.charAt(
                                                                                                                            Pt,
                                                                                                                        ),
                                                                                                                    )
                                                                                                                        ? ((s =
                                                                                                                              e.charAt(
                                                                                                                                  Pt,
                                                                                                                              )),
                                                                                                                          Pt++)
                                                                                                                        : ((s =
                                                                                                                              i),
                                                                                                                          0 ===
                                                                                                                              Bt &&
                                                                                                                              Ht(
                                                                                                                                  Q,
                                                                                                                              )));
                                                                                                            o !== i
                                                                                                                ? (r =
                                                                                                                      a =
                                                                                                                          [
                                                                                                                              a,
                                                                                                                              o,
                                                                                                                          ])
                                                                                                                : ((Pt =
                                                                                                                      r),
                                                                                                                  (r =
                                                                                                                      i));
                                                                                                        } else
                                                                                                            ((Pt = r),
                                                                                                                (r =
                                                                                                                    i));
                                                                                                        (r === i &&
                                                                                                            (kt.test(
                                                                                                                e.charAt(
                                                                                                                    Pt,
                                                                                                                ),
                                                                                                            )
                                                                                                                ? ((r =
                                                                                                                      e.charAt(
                                                                                                                          Pt,
                                                                                                                      )),
                                                                                                                  Pt++)
                                                                                                                : ((r =
                                                                                                                      i),
                                                                                                                  0 ===
                                                                                                                      Bt &&
                                                                                                                      Ht(
                                                                                                                          Rt,
                                                                                                                      ))),
                                                                                                            r !== i
                                                                                                                ? (t =
                                                                                                                      n =
                                                                                                                          [
                                                                                                                              n,
                                                                                                                              r,
                                                                                                                          ])
                                                                                                                : ((Pt =
                                                                                                                      t),
                                                                                                                  (t =
                                                                                                                      i)));
                                                                                                    } else
                                                                                                        ((Pt = t),
                                                                                                            (t = i));
                                                                                                    return (
                                                                                                        t !== i &&
                                                                                                            (t = Nt(t)),
                                                                                                        t
                                                                                                    );
                                                                                                })()),
                                                                                                c === i && (c = null),
                                                                                                c !== i
                                                                                                    ? (93 ===
                                                                                                      e.charCodeAt(Pt)
                                                                                                          ? ((u = A),
                                                                                                            Pt++)
                                                                                                          : ((u = i),
                                                                                                            0 === Bt &&
                                                                                                                Ht(_)),
                                                                                                      u !== i
                                                                                                          ? (t = n =
                                                                                                                [
                                                                                                                    n,
                                                                                                                    r,
                                                                                                                    a,
                                                                                                                    o,
                                                                                                                    s,
                                                                                                                    l,
                                                                                                                    c,
                                                                                                                    u,
                                                                                                                ])
                                                                                                          : ((Pt = t),
                                                                                                            (t = i)))
                                                                                                    : ((Pt = t),
                                                                                                      (t = i)))
                                                                                              : ((Pt = t), (t = i)))
                                                                                        : ((Pt = t), (t = i)))
                                                                                  : ((Pt = t), (t = i)))
                                                                            : ((Pt = t), (t = i)))
                                                                      : ((Pt = t), (t = i)))
                                                                : ((Pt = t), (t = i)),
                                                            t !== i && (t = C(t)),
                                                            t
                                                        );
                                                    })()),
                                                    t === i && (t = Yt())),
                                                t !== i && (t = g(t)),
                                                t
                                            );
                                        })()),
                                        n !== i)
                                    ) {
                                        for (r = [], a = Ut(); a !== i; ) (r.push(a), (a = Ut()));
                                        if (r !== i) {
                                            for (
                                                a = [],
                                                    o = Pt,
                                                    (l = Vt()) === i && (l = null),
                                                    l !== i && (c = Zt()) !== i
                                                        ? (o = l = [l, c])
                                                        : ((Pt = o), (o = i));
                                                o !== i;
                                            )
                                                (a.push(o),
                                                    (o = Pt),
                                                    (l = Vt()) === i && (l = null),
                                                    l !== i && (c = Zt()) !== i
                                                        ? (o = l = [l, c])
                                                        : ((Pt = o), (o = i)));
                                            if (a !== i) {
                                                for (o = [], l = Ut(); l !== i; ) (o.push(l), (l = Ut()));
                                                if (o !== i)
                                                    if (((l = Vt()) === i && (l = null), l !== i))
                                                        if (((c = jt()) === i && (c = null), c !== i)) {
                                                            for (u = [], h = Ut(); h !== i; ) (u.push(h), (h = Ut()));
                                                            u !== i
                                                                ? (t = n = [n, r, a, o, l, c, u])
                                                                : ((Pt = t), (t = i));
                                                        } else ((Pt = t), (t = i));
                                                    else ((Pt = t), (t = i));
                                                else ((Pt = t), (t = i));
                                            } else ((Pt = t), (t = i));
                                        } else ((Pt = t), (t = i));
                                    } else ((Pt = t), (t = i));
                                    return (t !== i && (t = s(t)), t);
                                }
                                function Ut() {
                                    var t, n, r, a, o;
                                    return (
                                        Pt,
                                        (t = Pt),
                                        40 === e.charCodeAt(Pt) ? ((n = l), Pt++) : ((n = i), 0 === Bt && Ht(c)),
                                        n !== i
                                            ? ((r = Vt()) === i && (r = null),
                                              r !== i && (a = jt()) !== i
                                                  ? (41 === e.charCodeAt(Pt)
                                                        ? ((o = u), Pt++)
                                                        : ((o = i), 0 === Bt && Ht(h)),
                                                    o !== i ? (t = n = [n, r, a, o]) : ((Pt = t), (t = i)))
                                                  : ((Pt = t), (t = i)))
                                            : ((Pt = t), (t = i)),
                                        t !== i && (t = d(t)),
                                        t
                                    );
                                }
                                function Vt() {
                                    var t;
                                    return (
                                        Pt,
                                        f.test(e.charAt(Pt))
                                            ? ((t = e.charAt(Pt)), Pt++)
                                            : ((t = i), 0 === Bt && Ht(p)),
                                        t !== i && (t = m(t)),
                                        t
                                    );
                                }
                                function Wt() {
                                    var t;
                                    return (
                                        Pt,
                                        B.test(e.charAt(Pt))
                                            ? ((t = e.charAt(Pt)), Pt++)
                                            : ((t = i), 0 === Bt && Ht(O)),
                                        t !== i && (t = g(t)),
                                        t
                                    );
                                }
                                function Yt() {
                                    var t;
                                    return (
                                        Pt,
                                        42 === e.charCodeAt(Pt) ? ((t = F), Pt++) : ((t = i), 0 === Bt && Ht(z)),
                                        t !== i && (t = q(t)),
                                        t
                                    );
                                }
                                function Zt() {
                                    var t, n, r, a;
                                    return (
                                        Pt,
                                        (t = Pt),
                                        37 === e.charCodeAt(Pt) ? ((n = Y), Pt++) : ((n = i), 0 === Bt && Ht(Z)),
                                        n !== i
                                            ? (X.test(e.charAt(Pt))
                                                  ? ((r = e.charAt(Pt)), Pt++)
                                                  : ((r = i), 0 === Bt && Ht(G)),
                                              r !== i
                                                  ? (K.test(e.charAt(Pt))
                                                        ? ((a = e.charAt(Pt)), Pt++)
                                                        : ((a = i), 0 === Bt && Ht(Q)),
                                                    a !== i ? (t = n = [n, r, a]) : ((Pt = t), (t = i)))
                                                  : ((Pt = t), (t = i)))
                                            : ((Pt = t), (t = i)),
                                        t === i &&
                                            (K.test(e.charAt(Pt))
                                                ? ((t = e.charAt(Pt)), Pt++)
                                                : ((t = i), 0 === Bt && Ht(Q))),
                                        t !== i && (t = J(t)),
                                        t
                                    );
                                }
                                if ((r = o()) !== i && Pt === e.length) return r;
                                throw (
                                    r !== i && Pt < e.length && Ht({ type: "end" }),
                                    (function (e, n, r) {
                                        return new t(t.buildMessage(e, n), e, n, r);
                                    })(
                                        Et,
                                        It < e.length ? e.charAt(It) : null,
                                        It < e.length ? qt(It, It + 1) : qt(It, It),
                                    )
                                );
                            },
                        }
                    );
                })();
            },
            421: (t, e, n) => {
                const r = n(348),
                    i = n(614),
                    a = (n(152), n(333));
                class o {
                    constructor(t) {
                        ((this.id = null),
                            (this.members = t),
                            (this.edges = []),
                            (this.insiders = []),
                            (this.neighbours = []),
                            (this.positioned = !1),
                            (this.center = new i(0, 0)),
                            (this.rings = []),
                            (this.isBridged = !1),
                            (this.isPartOfBridged = !1),
                            (this.isSpiro = !1),
                            (this.isFused = !1),
                            (this.centralAngle = 0),
                            (this.canFlip = !0));
                    }
                    clone() {
                        let t = new o(this.members);
                        return (
                            (t.id = this.id),
                            (t.insiders = r.clone(this.insiders)),
                            (t.neighbours = r.clone(this.neighbours)),
                            (t.positioned = this.positioned),
                            (t.center = this.center.clone()),
                            (t.rings = r.clone(this.rings)),
                            (t.isBridged = this.isBridged),
                            (t.isPartOfBridged = this.isPartOfBridged),
                            (t.isSpiro = this.isSpiro),
                            (t.isFused = this.isFused),
                            (t.centralAngle = this.centralAngle),
                            (t.canFlip = this.canFlip),
                            t
                        );
                    }
                    getSize() {
                        return this.members.length;
                    }
                    getPolygon(t) {
                        let e = [];
                        for (let n = 0; n < this.members.length; n++) e.push(t[this.members[n]].position);
                        return e;
                    }
                    getAngle() {
                        return Math.PI - this.centralAngle;
                    }
                    eachMember(t, e, n, r) {
                        let i = (n = n || 0 === n ? n : this.members[0]),
                            a = 0;
                        for (; null != i && a < 100; ) {
                            let o = i;
                            (e(o), (i = t[i].getNextInRing(t, this.id, r)), (r = o), i == n && (i = null), a++);
                        }
                    }
                    getOrderedNeighbours(t) {
                        let e = Array(this.neighbours.length);
                        for (let n = 0; n < this.neighbours.length; n++) {
                            let r = a.getVertices(t, this.id, this.neighbours[n]);
                            e[n] = { n: r.length, neighbour: this.neighbours[n] };
                        }
                        return (
                            e.sort(function (t, e) {
                                return e.n - t.n;
                            }),
                            e
                        );
                    }
                    isBenzeneLike(t) {
                        let e = this.getDoubleBondCount(t),
                            n = this.members.length;
                        return (3 === e && 6 === n) || (2 === e && 5 === n);
                    }
                    getDoubleBondCount(t) {
                        let e = 0;
                        for (let n = 0; n < this.members.length; n++) {
                            let r = t[this.members[n]].value;
                            ("=" !== r.bondType && "=" !== r.branchBond) || e++;
                        }
                        return e;
                    }
                    contains(t) {
                        for (let e = 0; e < this.members.length; e++) if (this.members[e] == t) return !0;
                        return !1;
                    }
                }
                t.exports = o;
            },
            333: (t, e, n) => {
                (n(152),
                    n(421),
                    (t.exports = class {
                        constructor(t, e) {
                            ((this.id = null),
                                (this.firstRingId = t.id),
                                (this.secondRingId = e.id),
                                (this.vertices = new Set()));
                            for (var n = 0; n < t.members.length; n++) {
                                let r = t.members[n];
                                for (let t = 0; t < e.members.length; t++) r === e.members[t] && this.addVertex(r);
                            }
                        }
                        addVertex(t) {
                            this.vertices.add(t);
                        }
                        updateOther(t, e) {
                            this.firstRingId === e ? (this.secondRingId = t) : (this.firstRingId = t);
                        }
                        containsRing(t) {
                            return this.firstRingId === t || this.secondRingId === t;
                        }
                        isBridge(t) {
                            if (this.vertices.size > 2) return !0;
                            for (let e of this.vertices) if (t[e].value.rings.length > 2) return !0;
                            return !1;
                        }
                        static isBridge(t, e, n, r) {
                            let i = null;
                            for (let a = 0; a < t.length; a++)
                                if (
                                    ((i = t[a]),
                                    (i.firstRingId === n && i.secondRingId === r) ||
                                        (i.firstRingId === r && i.secondRingId === n))
                                )
                                    return i.isBridge(e);
                            return !1;
                        }
                        static getNeighbours(t, e) {
                            let n = [];
                            for (let r = 0; r < t.length; r++) {
                                let i = t[r];
                                i.firstRingId === e
                                    ? n.push(i.secondRingId)
                                    : i.secondRingId === e && n.push(i.firstRingId);
                            }
                            return n;
                        }
                        static getVertices(t, e, n) {
                            for (let r = 0; r < t.length; r++) {
                                let i = t[r];
                                if (
                                    (i.firstRingId === e && i.secondRingId === n) ||
                                    (i.firstRingId === n && i.secondRingId === e)
                                )
                                    return [...i.vertices];
                            }
                        }
                    }));
            },
            688: (t, e, n) => {
                const r = n(707);
                class i {
                    static getRings(t) {
                        let e = t.getComponentsAdjacencyMatrix();
                        if (0 === e.length) return null;
                        let n = r.getConnectedComponents(e),
                            a = Array();
                        for (var o = 0; o < n.length; o++) {
                            let e = n[o],
                                r = t.getSubgraphAdjacencyMatrix([...e]),
                                c = new Uint16Array(r.length),
                                u = new Uint16Array(r.length);
                            for (var s = 0; s < r.length; s++) {
                                ((u[s] = 0), (c[s] = 0));
                                for (var l = 0; l < r[s].length; l++) c[s] += r[s][l];
                            }
                            let h = 0;
                            for (s = 0; s < r.length; s++) for (l = s + 1; l < r.length; l++) h += r[s][l];
                            let d = h - r.length + 1,
                                g = !0;
                            for (s = 0; s < c.length; s++) 3 !== c[s] && (g = !1);
                            if ((g && (d = 2 + h - r.length), 1 === d)) {
                                a.push([...e]);
                                continue;
                            }
                            let { d: f, pe: p, pe_prime: m } = i.getPathIncludedDistanceMatrices(r),
                                v = i.getRingCandidates(f, p, m),
                                y = i.getSSSR(v, f, r, p, m, c, u, d);
                            for (s = 0; s < y.length; s++) {
                                let t = Array(y[s].size),
                                    n = 0;
                                for (let r of y[s]) t[n++] = e[r];
                                a.push(t);
                            }
                        }
                        return a;
                    }
                    static matrixToString(t) {
                        let e = "";
                        for (var n = 0; n < t.length; n++) {
                            for (var r = 0; r < t[n].length; r++) e += t[n][r] + " ";
                            e += "\n";
                        }
                        return e;
                    }
                    static getPathIncludedDistanceMatrices(t) {
                        let e = t.length,
                            n = Array(e),
                            r = Array(e),
                            i = Array(e);
                        for (var a = 0, o = 0, s = 0, l = e; l--; ) {
                            ((n[l] = Array(e)), (r[l] = Array(e)), (i[l] = Array(e)));
                            for (var c = e; c--; )
                                ((n[l][c] = l === c || 1 === t[l][c] ? t[l][c] : Number.POSITIVE_INFINITY),
                                    1 === n[l][c] ? (r[l][c] = [[[l, c]]]) : (r[l][c] = Array()),
                                    (i[l][c] = Array()));
                        }
                        for (var u = e; u--; )
                            for (l = e; l--; )
                                for (c = e; c--; ) {
                                    const t = n[l][c],
                                        e = n[l][u] + n[u][c];
                                    if (t > e) {
                                        if (t === e + 1)
                                            for (i[l][c] = [r[l][c].length], a = r[l][c].length; a--; )
                                                for (i[l][c][a] = [r[l][c][a].length], o = r[l][c][a].length; o--; )
                                                    for (
                                                        i[l][c][a][o] = [r[l][c][a][o].length],
                                                            s = r[l][c][a][o].length;
                                                        s--;
                                                    )
                                                        i[l][c][a][o][s] = [r[l][c][a][o][0], r[l][c][a][o][1]];
                                        else i[l][c] = Array();
                                        for (n[l][c] = e, r[l][c] = [[]], a = r[l][u][0].length; a--; )
                                            r[l][c][0].push(r[l][u][0][a]);
                                        for (a = r[u][c][0].length; a--; ) r[l][c][0].push(r[u][c][0][a]);
                                    } else if (t === e) {
                                        if (r[l][u].length && r[u][c].length)
                                            if (r[l][c].length) {
                                                let t = Array();
                                                for (a = r[l][u][0].length; a--; ) t.push(r[l][u][0][a]);
                                                for (a = r[u][c][0].length; a--; ) t.push(r[u][c][0][a]);
                                                r[l][c].push(t);
                                            } else {
                                                let t = Array();
                                                for (a = r[l][u][0].length; a--; ) t.push(r[l][u][0][a]);
                                                for (a = r[u][c][0].length; a--; ) t.push(r[u][c][0][a]);
                                                r[l][c][0] = t;
                                            }
                                    } else if (t === e - 1)
                                        if (i[l][c].length) {
                                            let t = Array();
                                            for (a = r[l][u][0].length; a--; ) t.push(r[l][u][0][a]);
                                            for (a = r[u][c][0].length; a--; ) t.push(r[u][c][0][a]);
                                            i[l][c].push(t);
                                        } else {
                                            let t = Array();
                                            for (a = r[l][u][0].length; a--; ) t.push(r[l][u][0][a]);
                                            for (a = r[u][c][0].length; a--; ) t.push(r[u][c][0][a]);
                                            i[l][c][0] = t;
                                        }
                                }
                        return { d: n, pe: r, pe_prime: i };
                    }
                    static getRingCandidates(t, e, n) {
                        let r = t.length,
                            i = Array(),
                            a = 0;
                        for (let o = 0; o < r; o++)
                            for (let s = 0; s < r; s++)
                                0 === t[o][s] ||
                                    (1 === e[o][s].length && 0 === n[o][s]) ||
                                    ((a = 0 !== n[o][s].length ? 2 * (t[o][s] + 0.5) : 2 * t[o][s]),
                                    a !== 1 / 0 && i.push([a, e[o][s], n[o][s]]));
                        return (
                            i.sort(function (t, e) {
                                return t[0] - e[0];
                            }),
                            i
                        );
                    }
                    static getSSSR(t, e, n, r, a, o, s, l) {
                        let c = Array(),
                            u = Array();
                        for (let e = 0; e < t.length; e++)
                            if (t[e][0] % 2 != 0)
                                for (let r = 0; r < t[e][2].length; r++) {
                                    let a = t[e][1][0].concat(t[e][2][r]);
                                    for (var h = 0; h < a.length; h++)
                                        a[h][0].constructor === Array && (a[h] = a[h][0]);
                                    let d = i.bondsToAtoms(a);
                                    if (
                                        (i.getBondCount(d, n) !== d.size ||
                                            i.pathSetsContain(c, d, a, u, o, s) ||
                                            (c.push(d), (u = u.concat(a))),
                                        c.length > l)
                                    )
                                        return c;
                                }
                            else
                                for (let r = 0; r < t[e][1].length - 1; r++) {
                                    let a = t[e][1][r].concat(t[e][1][r + 1]);
                                    for (h = 0; h < a.length; h++) a[h][0].constructor === Array && (a[h] = a[h][0]);
                                    let d = i.bondsToAtoms(a);
                                    if (
                                        (i.getBondCount(d, n) !== d.size ||
                                            i.pathSetsContain(c, d, a, u, o, s) ||
                                            (c.push(d), (u = u.concat(a))),
                                        c.length > l)
                                    )
                                        return c;
                                }
                        return c;
                    }
                    static getEdgeCount(t) {
                        let e = 0,
                            n = t.length;
                        for (var r = n - 1; r--; ) for (var i = n; i--; ) 1 === t[r][i] && e++;
                        return e;
                    }
                    static getEdgeList(t) {
                        let e = t.length,
                            n = Array();
                        for (var r = e - 1; r--; ) for (var i = e; i--; ) 1 === t[r][i] && n.push([r, i]);
                        return n;
                    }
                    static bondsToAtoms(t) {
                        let e = new Set();
                        for (var n = t.length; n--; ) (e.add(t[n][0]), e.add(t[n][1]));
                        return e;
                    }
                    static getBondCount(t, e) {
                        let n = 0;
                        for (let r of t) for (let i of t) r !== i && (n += e[r][i]);
                        return n / 2;
                    }
                    static pathSetsContain(t, e, n, r, a, o) {
                        for (var s = t.length; s--; ) {
                            if (i.isSupersetOf(e, t[s])) return !0;
                            if (t[s].size === e.size && i.areSetsEqual(t[s], e)) return !0;
                        }
                        let l = 0,
                            c = !1;
                        for (s = n.length; s--; )
                            for (var u = r.length; u--; )
                                (((n[s][0] === r[u][0] && n[s][1] === r[u][1]) ||
                                    (n[s][1] === r[u][0] && n[s][0] === r[u][1])) &&
                                    l++,
                                    l === n.length && (c = !0));
                        let h = !1;
                        if (c)
                            for (let t of e)
                                if (o[t] < a[t]) {
                                    h = !0;
                                    break;
                                }
                        if (c && !h) return !0;
                        for (let t of e) o[t]++;
                        return !1;
                    }
                    static areSetsEqual(t, e) {
                        if (t.size !== e.size) return !1;
                        for (let n of t) if (!e.has(n)) return !1;
                        return !0;
                    }
                    static isSupersetOf(t, e) {
                        for (var n of e) if (!t.has(n)) return !1;
                        return !0;
                    }
                }
                t.exports = i;
            },
            614: (t) => {
                class e {
                    constructor(t, e) {
                        0 == arguments.length
                            ? ((this.x = 0), (this.y = 0))
                            : 1 == arguments.length
                              ? ((this.x = t.x), (this.y = t.y))
                              : ((this.x = t), (this.y = e));
                    }
                    clone() {
                        return new e(this.x, this.y);
                    }
                    toString() {
                        return "(" + this.x + "," + this.y + ")";
                    }
                    add(t) {
                        return ((this.x += t.x), (this.y += t.y), this);
                    }
                    subtract(t) {
                        return ((this.x -= t.x), (this.y -= t.y), this);
                    }
                    divide(t) {
                        return ((this.x /= t), (this.y /= t), this);
                    }
                    multiply(t) {
                        return ((this.x *= t.x), (this.y *= t.y), this);
                    }
                    multiplyScalar(t) {
                        return ((this.x *= t), (this.y *= t), this);
                    }
                    invert() {
                        return ((this.x = -this.x), (this.y = -this.y), this);
                    }
                    angle() {
                        return Math.atan2(this.y, this.x);
                    }
                    distance(t) {
                        return Math.sqrt((t.x - this.x) * (t.x - this.x) + (t.y - this.y) * (t.y - this.y));
                    }
                    distanceSq(t) {
                        return (t.x - this.x) * (t.x - this.x) + (t.y - this.y) * (t.y - this.y);
                    }
                    clockwise(t) {
                        let e = this.y * t.x,
                            n = this.x * t.y;
                        return e > n ? -1 : e === n ? 0 : 1;
                    }
                    relativeClockwise(t, e) {
                        let n = (this.y - t.y) * (e.x - t.x),
                            r = (this.x - t.x) * (e.y - t.y);
                        return n > r ? -1 : n === r ? 0 : 1;
                    }
                    rotate(t) {
                        let n = new e(0, 0),
                            r = Math.cos(t),
                            i = Math.sin(t);
                        return (
                            (n.x = this.x * r - this.y * i),
                            (n.y = this.x * i + this.y * r),
                            (this.x = n.x),
                            (this.y = n.y),
                            this
                        );
                    }
                    rotateAround(t, e) {
                        let n = Math.sin(t),
                            r = Math.cos(t);
                        ((this.x -= e.x), (this.y -= e.y));
                        let i = this.x * r - this.y * n,
                            a = this.x * n + this.y * r;
                        return ((this.x = i + e.x), (this.y = a + e.y), this);
                    }
                    rotateTo(t, n, r = 0) {
                        ((this.x += 0.001), (this.y -= 0.001));
                        let i = e.subtract(this, n),
                            a = e.subtract(t, n),
                            o = e.angle(a, i);
                        return (this.rotateAround(o + r, n), this);
                    }
                    rotateAwayFrom(t, e, n) {
                        this.rotateAround(n, e);
                        let r = this.distanceSq(t);
                        (this.rotateAround(-2 * n, e), this.distanceSq(t) < r && this.rotateAround(2 * n, e));
                    }
                    getRotateAwayFromAngle(t, e, n) {
                        let r = this.clone();
                        r.rotateAround(n, e);
                        let i = r.distanceSq(t);
                        return (r.rotateAround(-2 * n, e), r.distanceSq(t) < i ? n : -n);
                    }
                    getRotateTowardsAngle(t, e, n) {
                        let r = this.clone();
                        r.rotateAround(n, e);
                        let i = r.distanceSq(t);
                        return (r.rotateAround(-2 * n, e), r.distanceSq(t) > i ? n : -n);
                    }
                    getRotateToAngle(t, n) {
                        let r = e.subtract(this, n),
                            i = e.subtract(t, n),
                            a = e.angle(i, r);
                        return Number.isNaN(a) ? 0 : a;
                    }
                    isInPolygon(t) {
                        let e = !1;
                        for (let n = 0, r = t.length - 1; n < t.length; r = n++)
                            t[n].y > this.y != t[r].y > this.y &&
                                this.x < ((t[r].x - t[n].x) * (this.y - t[n].y)) / (t[r].y - t[n].y) + t[n].x &&
                                (e = !e);
                        return e;
                    }
                    length() {
                        return Math.sqrt(this.x * this.x + this.y * this.y);
                    }
                    lengthSq() {
                        return this.x * this.x + this.y * this.y;
                    }
                    normalize() {
                        return (this.divide(this.length()), this);
                    }
                    normalized() {
                        return e.divideScalar(this, this.length());
                    }
                    whichSide(t, e) {
                        return (this.x - t.x) * (e.y - t.y) - (this.y - t.y) * (e.x - t.x);
                    }
                    sameSideAs(t, e, n) {
                        let r = this.whichSide(t, e),
                            i = n.whichSide(t, e);
                        return (r < 0 && i < 0) || (0 == r && 0 == i) || (r > 0 && i > 0);
                    }
                    static add(t, n) {
                        return new e(t.x + n.x, t.y + n.y);
                    }
                    static subtract(t, n) {
                        return new e(t.x - n.x, t.y - n.y);
                    }
                    static multiply(t, n) {
                        return new e(t.x * n.x, t.y * n.y);
                    }
                    static multiplyScalar(t, n) {
                        return new e(t.x, t.y).multiplyScalar(n);
                    }
                    static midpoint(t, n) {
                        return new e((t.x + n.x) / 2, (t.y + n.y) / 2);
                    }
                    static normals(t, n) {
                        let r = e.subtract(n, t);
                        return [new e(-r.y, r.x), new e(r.y, -r.x)];
                    }
                    static units(t, n) {
                        let r = e.subtract(n, t);
                        return [new e(-r.y, r.x).normalize(), new e(r.y, -r.x).normalize()];
                    }
                    static divide(t, n) {
                        return new e(t.x / n.x, t.y / n.y);
                    }
                    static divideScalar(t, n) {
                        return new e(t.x / n, t.y / n);
                    }
                    static dot(t, e) {
                        return t.x * e.x + t.y * e.y;
                    }
                    static angle(t, n) {
                        let r = e.dot(t, n);
                        return Math.acos(r / (t.length() * n.length()));
                    }
                    static threePointangle(t, n, r) {
                        let i = e.subtract(n, t),
                            a = e.subtract(r, n),
                            o = t.distance(n),
                            s = n.distance(r);
                        return Math.acos(e.dot(i, a) / (o * s));
                    }
                    static scalarProjection(t, n) {
                        let r = n.normalized();
                        return e.dot(t, r);
                    }
                    static averageDirection(t) {
                        let n = new e(0, 0);
                        for (var r = 0; r < t.length; r++) {
                            let e = t[r];
                            n.add(e);
                        }
                        return n.normalize();
                    }
                }
                t.exports = e;
            },
            152: (t, e, n) => {
                const r = n(474),
                    i = n(348),
                    a = n(614);
                n(427);
                class o {
                    constructor(t, e = 0, n = 0) {
                        ((this.id = null),
                            (this.value = t),
                            (this.position = new a(e || 0, n || 0)),
                            (this.previousPosition = new a(0, 0)),
                            (this.parentVertexId = null),
                            (this.children = Array()),
                            (this.spanningTreeChildren = Array()),
                            (this.edges = Array()),
                            (this.positioned = !1),
                            (this.angle = null),
                            (this.dir = 1),
                            (this.neighbourCount = 0),
                            (this.neighbours = Array()),
                            (this.neighbouringElements = Array()),
                            (this.forcePositioned = !1));
                    }
                    setPosition(t, e) {
                        ((this.position.x = t), (this.position.y = e));
                    }
                    setPositionFromVector(t) {
                        ((this.position.x = t.x), (this.position.y = t.y));
                    }
                    addChild(t) {
                        (this.children.push(t), this.neighbours.push(t), this.neighbourCount++);
                    }
                    addRingbondChild(t, e) {
                        if ((this.children.push(t), this.value.bracket)) {
                            let n = 1;
                            (0 === this.id && 0 === this.value.bracket.hcount && (n = 0),
                                1 === this.value.bracket.hcount && 0 === e && (n = 2),
                                1 === this.value.bracket.hcount && 1 === e && (n = this.neighbours.length < 3 ? 2 : 3),
                                null === this.value.bracket.hcount && 0 === e && (n = 1),
                                null === this.value.bracket.hcount &&
                                    1 === e &&
                                    (n = this.neighbours.length < 3 ? 1 : 2),
                                this.neighbours.splice(n, 0, t));
                        } else this.neighbours.push(t);
                        this.neighbourCount++;
                    }
                    setParentVertexId(t) {
                        (this.neighbourCount++, (this.parentVertexId = t), this.neighbours.push(t));
                    }
                    isTerminal() {
                        return (
                            !!this.value.hasAttachedPseudoElements ||
                            (null === this.parentVertexId && this.children.length < 2) ||
                            0 === this.children.length
                        );
                    }
                    clone() {
                        let t = new o(this.value, this.position.x, this.position.y);
                        return (
                            (t.id = this.id),
                            (t.previousPosition = new a(this.previousPosition.x, this.previousPosition.y)),
                            (t.parentVertexId = this.parentVertexId),
                            (t.children = i.clone(this.children)),
                            (t.spanningTreeChildren = i.clone(this.spanningTreeChildren)),
                            (t.edges = i.clone(this.edges)),
                            (t.positioned = this.positioned),
                            (t.angle = this.angle),
                            (t.forcePositioned = this.forcePositioned),
                            t
                        );
                    }
                    equals(t) {
                        return this.id === t.id;
                    }
                    getAngle(t = null, e = !1) {
                        let n = null;
                        return (
                            (n = t ? a.subtract(this.position, t) : a.subtract(this.position, this.previousPosition)),
                            e ? r.toDeg(n.angle()) : n.angle()
                        );
                    }
                    getTextDirection(t) {
                        let e = this.getDrawnNeighbours(t),
                            n = Array();
                        for (let r = 0; r < e.length; r++) n.push(this.getAngle(t[e[r]].position));
                        let i = r.meanAngle(n),
                            a = Math.PI / 2;
                        return (
                            (i = Math.round(Math.round(i / a) * a)),
                            2 === i
                                ? "down"
                                : -2 === i
                                  ? "up"
                                  : 0 === i || -0 === i
                                    ? "right"
                                    : 3 === i || -3 === i
                                      ? "left"
                                      : "down"
                        );
                    }
                    getNeighbours(t = null) {
                        if (null === t) return this.neighbours.slice();
                        let e = Array();
                        for (let n = 0; n < this.neighbours.length; n++)
                            this.neighbours[n] !== t && e.push(this.neighbours[n]);
                        return e;
                    }
                    getDrawnNeighbours(t) {
                        let e = Array();
                        for (let n = 0; n < this.neighbours.length; n++)
                            t[this.neighbours[n]].value.isDrawn && e.push(this.neighbours[n]);
                        return e;
                    }
                    getNeighbourCount() {
                        return this.neighbourCount;
                    }
                    getSpanningTreeNeighbours(t = null) {
                        let e = Array();
                        for (let n = 0; n < this.spanningTreeChildren.length; n++)
                            (void 0 !== t && t == this.spanningTreeChildren[n]) || e.push(this.spanningTreeChildren[n]);
                        return (
                            null != this.parentVertexId &&
                                ((void 0 !== t && t == this.parentVertexId) || e.push(this.parentVertexId)),
                            e
                        );
                    }
                    getNextInRing(t, e, n) {
                        let r = this.getNeighbours();
                        for (let a = 0; a < r.length; a++)
                            if (i.contains(t[r[a]].value.rings, { value: e }) && r[a] != n) return r[a];
                        return null;
                    }
                }
                t.exports = o;
            },
        },
        e = {};
    function n(r) {
        var i = e[r];
        if (void 0 !== i) return i.exports;
        var a = (e[r] = { exports: {} });
        return (t[r].call(a.exports, a, a.exports, n), a.exports);
    }
    ((n.d = (t, e) => {
        for (var r in e) n.o(e, r) && !n.o(t, r) && Object.defineProperty(t, r, { enumerable: !0, get: e[r] });
    }),
        (n.o = (t, e) => Object.prototype.hasOwnProperty.call(t, e)),
        (n.r = (t) => {
            ("undefined" != typeof Symbol &&
                Symbol.toStringTag &&
                Object.defineProperty(t, Symbol.toStringTag, { value: "Module" }),
                Object.defineProperty(t, "__esModule", { value: !0 }));
        }));
    var r = {};
    ((() => {
        "use strict";
        var t = r;
        (Object.defineProperty(t, "__esModule", { value: !0 }), (t.start = t.getAnchor = t.downloadSvg = void 0));
        const e = n(220),
            i = n(225),
            a = n(26),
            o = n(322),
            s = n(128),
            l = n(501),
            c = n(415),
            u = n(82),
            h = n(457),
            d = n(878),
            g = n(681),
            f = n(625),
            p = n(24),
            m = n(3);
        var v = n(415);
        Object.defineProperty(t, "downloadSvg", {
            enumerable: !0,
            get: function () {
                return v.downloadSvg;
            },
        });
        const y = "antismash.outputs.html.visualisers";
        let b = null,
            x = null,
            w = null;
        function A(t) {
            (t.preventDefault(), $("#downloadmenu").fadeToggle("fast", "linear"));
        }
        function _() {
            return window.location.hash.substring(1) || "overview";
        }
        function C() {
            setTimeout(() => {
                ($(".page").hide(), $(".empty-on-leave").empty(), $(".regbutton").removeClass("active"));
                const t = _();
                ($(`#${t}`).show(),
                    "overview" !== t &&
                        ($(`.regbutton.${t}`).addClass("active"),
                        void 0 !== b[t] && (0, m.drawRegion)(`${t}-svg`, b[t], 20),
                        $(`#${t}-details-svg`).length > 0 && (0, d.drawDomains)(t, x.nrpspks[t], 25),
                        $(`#${t} .clusterblast-selector`).change(),
                        t in w &&
                            ("antismash.modules.cluster_compare" in w[t] &&
                                (0, o.setComparisonData)(t, w[t]["antismash.modules.cluster_compare"], b[t]),
                            `${y}.bubble_view` in w[t] && (0, l.drawDomainBubbleData)(t, w[t][`${y}.bubble_view`]),
                            `${y}.generic_domains` in w[t] &&
                                (0, u.drawGenericDomains)(t, w[t][`${y}.generic_domains`], 25),
                            `${y}.gene_table` in w[t] && (0, h.initGeneTableHandler)(b[t], w[t][`${y}.gene_table`]),
                            "antismash.modules.tfbs_finder" in w[t] &&
                                (0, p.drawBindingSites)(t, w[t]["antismash.modules.tfbs_finder"])),
                        $(`#${t} .comparison-selector`).change(),
                        $(`#${t} * .body-details-header-active`).first().click()));
            }, 1);
        }
        function S() {
            const t = b.order,
                e = _();
            let n = "overview";
            if ("overview" === e) n = t[0];
            else {
                const r = t.indexOf(e);
                r !== t.length - 1 && (n = t[r + 1]);
            }
            ((window.location.href = `#${n}`), C());
        }
        function M() {
            const t = b.order,
                e = _();
            let n = "";
            if ("overview" === e) n = t[t.length - 1];
            else {
                const r = t.indexOf(e);
                0 !== r && (n = t[r - 1]);
            }
            ((window.location.href = `#${n}`), C());
        }
        function T(t, e) {
            const n = t + "-active",
                r = _(),
                i = $(`#${r}`).find(`.${t}`),
                a = $(`#${r}`).find(`.${n}`);
            i &&
                a &&
                a &&
                (e > 0
                    ? a.is(i.last())
                        ? i.first().trigger("click")
                        : a.next().trigger("click")
                    : a.is(i.first())
                      ? i.last().trigger("click")
                      : a.prev().trigger("click"));
        }
        function k(t) {
            const e = t.keyCode;
            (37 === e || 87 === e ? M() : (39 !== e && 69 !== e) || S(),
                "overview" !== _() &&
                    (65 === e
                        ? T("body-details-header", -1)
                        : 83 === e
                          ? T("body-details-header", 1)
                          : 68 === e
                            ? T("sidepanel-details-header", -1)
                            : 70 === e
                              ? T("sidepanel-details-header", 1)
                              : 77 === e && $("input.show-module-domains").first().click()));
        }
        function R(t) {
            t.preventDefault();
            const e = ($(this).attr("id") || "").replace(/-header/, "");
            if (!e) return;
            const n = $(`#${e}`);
            ("none" === n.css("display")
                ? $(this).text("Hide pHMM detection rules used")
                : $(this).text("Show pHMM detection rules used"),
                n.fadeToggle("fast", "linear"));
        }
        ((t.getAnchor = _),
            (t.start = function (t, n, r, o) {
                ((0, g.createRecordOverviews)(o),
                    (b = t),
                    (x = n),
                    (w = r),
                    document.addEventListener("keyup", k, !1),
                    $("#download").click(A),
                    $("#next-region").click(S),
                    $("#prev-region").click(M),
                    $(".regbutton")
                        .click(function () {
                            const t = $(this).children().first().attr("href");
                            void 0 !== t && ((window.location.href = t), C());
                        })
                        .mouseover(function () {
                            const t = ($(this).attr("class") || "").split(" ");
                            if (t.length < 2) return;
                            if ("separator" === t[1]) return;
                            const e = (function (t) {
                                    switch (t) {
                                        case "nrps":
                                            return "NRPS";
                                        case "t1pks":
                                            return "Type I PKS";
                                        case "t2pks":
                                            return "Type II PKS";
                                        case "t3pks":
                                            return "Type III PKS";
                                        case "t4pks":
                                            return "Type IV PKS";
                                        default:
                                            return t;
                                    }
                                })(t[1]),
                                n = $("#region-type");
                            (n.data("orig_text", n.text()), n.text(`${e}:`));
                        })
                        .mouseout(() => {
                            const t = $("#region-type");
                            t.text(t.data("orig_text"));
                        }),
                    $(".clusterblast-selector").change(function () {
                        const t = ($(this).attr("id") || "nonexistant").replace("-select", ""),
                            e = "" + $(this).val();
                        (e &&
                            $.get(
                                e,
                                (e) => {
                                    ($(`#${t}-svg`).html(e), (0, i.init)(`${t}-svg`));
                                },
                                "html",
                            ),
                            $(`#${t}-download`).off("click"),
                            $(`#${t}-download`).click(() => window.open("" + $(`#${t}-select`).val(), "_blank")));
                    }),
                    $(".comparison-selector").change(function () {
                        const t = ($(this).attr("id") || "nonexistant").replace("-selector", "");
                        $(`#${t}`).siblings().removeClass("comparison-container-active");
                        const e = $(this).attr("data-tag");
                        $(`#${t}-${$(this).val()}`)
                            .addClass("comparison-container-active")
                            .find(`.heat-row-${e}`)
                            .first()
                            .click();
                    }),
                    $(".cluster-rules-header").click(R),
                    C(),
                    $("input.overview-switch-compact")
                        .change(function () {
                            $(this).prop("checked")
                                ? ($("#single-record-tables").hide(), $("#compact-record-table").show())
                                : ($("#single-record-tables").show(), $("#compact-record-table").hide());
                        })
                        .trigger("change"),
                    $(".linked-row")
                        .off("click")
                        .click(function () {
                            const t = $(this).attr("data-anchor");
                            t && ((window.location.href = t), C());
                        }),
                    $("input.domains-selected-only").change(function () {
                        ($("input.domains-selected-only").prop("checked", $(this).prop("checked")),
                            (0, d.redrawDomains)(),
                            (0, u.redrawGenericDomains)());
                    }),
                    $("input.domains-expand-full").change(function () {
                        ($("input.domains-expand-full").prop("checked", $(this).prop("checked")),
                            (0, u.redrawGenericDomains)({ reset: !0, anchor: _() }));
                    }),
                    (0, d.createModuleHandlers)(),
                    (0, f.drawStructures)(),
                    (0, s.setupDetails)(t.order),
                    $(".help-icon")
                        .off("click")
                        .click(function () {
                            ($(this).toggleClass("active"), $(this).next().toggle());
                        }),
                    $(".help-tooltip")
                        .off("click")
                        .click(function () {
                            ($(this).hide(), $(this).prev().removeClass("active"));
                        }),
                    $(".clipboard-copy").off("click").click(e.copyToClipboard),
                    $(".collapser").click(a.toggleCollapserHandler),
                    (0, c.initDownloadButtons)());
            }));
    })(),
        (viewer = r));
})();
