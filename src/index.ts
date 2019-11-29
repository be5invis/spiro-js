import { DefaultBezierContext, IBezierContext, KnotCallback, Point } from "./context";
import { integrate_spiro, SpiroArc, SpiroK } from "./spiro-arc";

export * from "./context";
export { SpiroArc, SpiroK } from "./spiro-arc";

const MAX_DEPTH = 7;

export type PointType = "corner" | "open" | "open_end" | "left" | "right" | "g2" | "g4";
export type Knot<C> = Point & { type: PointType; af?: KnotCallback<C> };

class SpiroSeg<C> {
	x = 0;
	y = 0;
	type: PointType = "corner";
	bend_th = 0;
	ks: SpiroK = [0, 0, 0, 0];
	seg_ch = 0;
	seg_th = 0;
	l = 0;
	af?: KnotCallback<C>;
}

class BandMat {
	a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // band-diagonal matrix
	al = [0, 0, 0, 0, 0]; // lower part of band-diagonal decomposition
}

function copyArray<T>(from: T[], fromI: number, to: T[], toI: number, nElem: number) {
	for (let k = 0; k < nElem; k++) {
		to[toI + k] = from[fromI + k];
	}
}

function copyBandMatArray(
	from: BandMat[],
	fromI: number,
	to: BandMat[],
	toI: number,
	nElem: number
) {
	for (let i = 0; i < nElem; ++i) {
		to[i + toI].a = [...from[i + fromI].a];
		to[i + toI].al = [...from[i + fromI].al];
	}
}

function compute_ends(ks: SpiroK, ends: number[][], seg_ch: number) {
	const xy = integrate_spiro(ks);
	const ch = Math.hypot(xy[0], xy[1]);
	const th = Math.atan2(xy[1], xy[0]);
	const l = ch / seg_ch;

	const th_even = 0.5 * ks[0] + (1 / 48) * ks[2];
	const th_odd = 0.125 * ks[1] + (1 / 384) * ks[3] - th;
	ends[0][0] = th_even - th_odd;
	ends[1][0] = th_even + th_odd;

	const k0_even = l * (ks[0] + 0.125 * ks[2]);
	const k0_odd = l * (0.5 * ks[1] + (1 / 48) * ks[3]);
	ends[0][1] = k0_even - k0_odd;
	ends[1][1] = k0_even + k0_odd;

	const l2 = l * l;
	const k1_even = l2 * (ks[1] + 0.125 * ks[3]);
	const k1_odd = l2 * 0.5 * ks[2];
	ends[0][2] = k1_even - k1_odd;
	ends[1][2] = k1_even + k1_odd;

	const l3 = l2 * l;
	const k2_even = l3 * ks[2];
	const k2_odd = l3 * 0.5 * ks[3];
	ends[0][3] = k2_even - k2_odd;
	ends[1][3] = k2_even + k2_odd;

	return l;
}

function compute_pderivs<C>(s: SpiroSeg<C>, ends: number[][], derivs: number[][][], jinc: number) {
	const recip_d = 2e6;
	const delta = 1 / recip_d;
	let try_ks: SpiroK = [0, 0, 0, 0];
	let try_ends = [
		[0, 0, 0, 0],
		[0, 0, 0, 0]
	];

	compute_ends(s.ks, ends, s.seg_ch);
	for (let i = 0; i < jinc; i++) {
		for (let j = 0; j < 4; j++) try_ks[j] = s.ks[j];

		try_ks[i] += delta;
		compute_ends(try_ks, try_ends, s.seg_ch);
		for (let k = 0; k < 2; k++) {
			for (let j = 0; j < 4; j++) {
				derivs[j][k][i] = recip_d * (try_ends[k][j] - ends[k][j]);
			}
		}
	}
}

function mod_2pi(th: number) {
	const u = th / (2 * Math.PI);
	return 2 * Math.PI * (u - Math.floor(u + 0.5));
}

function setup_path<C>(src: Knot<C>[], n: number) {
	let n_seg = src[0].type === "open" ? n - 1 : n;
	let r: SpiroSeg<C>[] = [];
	for (let j = 0; j < n_seg + 1; j++) {
		r[j] = new SpiroSeg();
	}

	for (let i = 0; i < n_seg; i++) {
		r[i] = new SpiroSeg();
		r[i].x = src[i].x;
		r[i].y = src[i].y;
		r[i].type = src[i].type;
		r[i].af = src[i].af;
		r[i].ks = [0, 0, 0, 0];
	}
	r[n_seg] = new SpiroSeg();
	r[n_seg].x = src[n_seg % n].x;
	r[n_seg].y = src[n_seg % n].y;
	r[n_seg].type = src[n_seg % n].type;
	r[n_seg].af = src[n_seg % n].af;

	for (let i = 0; i < n_seg; i++) {
		let dx = r[i + 1].x - r[i].x;
		let dy = r[i + 1].y - r[i].y;
		r[i].seg_ch = Math.hypot(dx, dy);
		r[i].seg_th = Math.atan2(dy, dx);
	}

	let iLast = n_seg - 1;
	for (let i = 0; i < n_seg; i++) {
		if (r[i].type === "open" || r[i].type === "open_end" || r[i].type === "corner") {
			r[i].bend_th = 0;
		} else {
			r[i].bend_th = mod_2pi(r[i].seg_th - r[iLast].seg_th);
		}
		iLast = i;
	}
	return r;
}

function bandec11(m: BandMat[], perm: number[], n: number) {
	let i: number, j: number, k: number;
	let l: number;

	/* pack top triangle to the LEFT. */
	for (i = 0; i < 5; i++) {
		for (j = 0; j < i + 6; j++) {
			m[i].a[j] = m[i].a[j + 5 - i];
		}
		for (; j < 11; j++) m[i].a[j] = 0;
	}
	l = 5;
	for (k = 0; k < n; k++) {
		let pivot = k;
		let pivot_val = m[k].a[0];

		l = l < n ? l + 1 : n;

		for (j = k + 1; j < l; j++) {
			if (Math.abs(m[j].a[0]) > Math.abs(pivot_val)) {
				pivot_val = m[j].a[0];
				pivot = j;
			}
		}

		perm[k] = pivot;
		if (pivot !== k) {
			for (j = 0; j < 11; j++) {
				let tmp = m[k].a[j];
				m[k].a[j] = m[pivot].a[j];
				m[pivot].a[j] = tmp;
			}
		}

		if (Math.abs(pivot_val) < 1e-12) pivot_val = 1e-12;
		const pivot_scale = 1 / pivot_val;
		for (i = k + 1; i < l; i++) {
			let x = m[i].a[0] * pivot_scale;
			m[k].al[i - k - 1] = x;
			for (j = 1; j < 11; j++) {
				m[i].a[j - 1] = m[i].a[j] - x * m[k].a[j];
			}
			m[i].a[10] = 0;
		}
	}
}

function banbks11(m: BandMat[], perm: number[], v: number[], n: number) {
	let i: number, k: number, l: number;

	/* forward substitution */
	l = 5;
	for (k = 0; k < n; k++) {
		i = perm[k];
		if (i !== k) {
			let tmp = v[k];
			v[k] = v[i];
			v[i] = tmp;
		}
		if (l < n) l++;
		for (i = k + 1; i < l; i++) {
			v[i] -= m[k].al[i - k - 1] * v[k];
		}
	}

	/* back substitution */
	l = 1;
	for (i = n - 1; i >= 0; i--) {
		let x = v[i];
		for (k = 1; k < l; k++) {
			x -= m[i].a[k] * v[k + i];
		}
		v[i] = x / m[i].a[0];
		if (l < 11) l++;
	}
}

function compute_jinc(ty0: PointType, ty1: PointType) {
	if (ty0 === "g4" || ty1 === "g4" || ty0 === "right" || ty1 === "left") {
		return 4;
	} else if (ty0 === "g2" && ty1 === "g2") {
		return 2;
	} else if (
		((ty0 === "open" || ty0 === "corner" || ty0 === "left") && ty1 === "g2") ||
		(ty0 === "g2" && (ty1 === "open_end" || ty1 === "corner" || ty1 === "right"))
	) {
		return 1;
	} else {
		return 0;
	}
}

function count_vec<C>(s: SpiroSeg<C>[], nseg: number) {
	let n = 0;
	for (let i = 0; i < nseg; i++) {
		n += compute_jinc(s[i].type, s[i + 1].type);
	}
	return n;
}

function add_mat_line(
	m: BandMat[],
	v: number[],
	derivs: number[],
	x: number,
	y: number,
	j: number,
	jj: number,
	jinc: number,
	nmat: number
) {
	if (jj >= 0) {
		let joff = (j + 5 - jj + nmat) % nmat;
		if (nmat < 6) {
			joff = j + 5 - jj;
		} else if (nmat === 6) {
			joff = 2 + ((j + 3 - jj + nmat) % nmat);
		}
		v[jj] += x;
		for (let k = 0; k < jinc; k++) {
			m[jj].a[joff + k] += y * derivs[k];
		}
	}
}

function spiro_iter<C>(s: SpiroSeg<C>[], m: BandMat[], perm: number[], v: number[], n: number) {
	const cyclic = s[0].type !== "open" && s[0].type !== "corner";
	let i: number, j: number, jj: number;
	let nmat = count_vec(s, n);
	let norm: number;
	let n_invert: number;

	for (i = 0; i < nmat; i++) {
		v[i] = 0;
		m[i] = new BandMat();
	}

	j = 0;
	if (s[0].type === "g4") {
		jj = nmat - 2;
	} else if (s[0].type === "g2") {
		jj = nmat - 1;
	} else {
		jj = 0;
	}
	for (i = 0; i < n; i++) {
		const ty0 = s[i].type;
		const ty1 = s[i + 1].type;
		const jinc = compute_jinc(ty0, ty1);
		const th = s[i].bend_th;
		let ends = [
			[0, 0, 0, 0],
			[0, 0, 0, 0]
		];
		let derivs = [
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			]
		];
		let jthl = -1,
			jk0l = -1,
			jk1l = -1,
			jk2l = -1;
		let jthr = -1,
			jk0r = -1,
			jk1r = -1,
			jk2r = -1;

		compute_pderivs(s[i], ends, derivs, jinc);

		/* constraints crossing LEFT */
		if (ty0 === "g4" || ty0 === "g2" || ty0 === "left" || ty0 === "right") {
			jthl = jj++;
			jj %= nmat;
			jk0l = jj++;
		}
		if (ty0 === "g4") {
			jj %= nmat;
			jk1l = jj++;
			jk2l = jj++;
		}

		/* constraints on LEFT */
		if ((ty0 === "left" || ty0 === "corner" || ty0 === "open" || ty0 === "g2") && jinc === 4) {
			if (ty0 !== "g2") jk1l = jj++;
			jk2l = jj++;
		}

		/* constraints on RIGHT */
		if (
			(ty1 === "right" || ty1 === "corner" || ty1 === "open_end" || ty1 === "g2") &&
			jinc === 4
		) {
			if (ty1 !== "g2") jk1r = jj++;
			jk2r = jj++;
		}

		/* constraints crossing RIGHT */
		if (ty1 === "g4" || ty1 === "g2" || ty1 === "left" || ty1 === "right") {
			jthr = jj;
			jk0r = (jj + 1) % nmat;
		}
		if (ty1 === "g4") {
			jk1r = (jj + 2) % nmat;
			jk2r = (jj + 3) % nmat;
		}

		add_mat_line(m, v, derivs[0][0], th - ends[0][0], 1, j, jthl, jinc, nmat);
		add_mat_line(m, v, derivs[1][0], ends[0][1], -1, j, jk0l, jinc, nmat);
		add_mat_line(m, v, derivs[2][0], ends[0][2], -1, j, jk1l, jinc, nmat);
		add_mat_line(m, v, derivs[3][0], ends[0][3], -1, j, jk2l, jinc, nmat);
		add_mat_line(m, v, derivs[0][1], -ends[1][0], 1, j, jthr, jinc, nmat);
		add_mat_line(m, v, derivs[1][1], -ends[1][1], 1, j, jk0r, jinc, nmat);
		add_mat_line(m, v, derivs[2][1], -ends[1][2], 1, j, jk1r, jinc, nmat);
		add_mat_line(m, v, derivs[3][1], -ends[1][3], 1, j, jk2r, jinc, nmat);
		j += jinc;
	}
	if (cyclic) {
		copyBandMatArray(m, 0, m, nmat, nmat);
		copyBandMatArray(m, 0, m, 2 * nmat, nmat);
		copyArray(v, 0, v, nmat, nmat);
		copyArray(v, 0, v, 2 * nmat, nmat);
		n_invert = 3 * nmat;
		j = nmat;
	} else {
		n_invert = nmat;
		j = 0;
	}
	bandec11(m, perm, n_invert);
	banbks11(m, perm, v, n_invert);
	norm = 0;
	for (i = 0; i < n; i++) {
		const ty0 = s[i].type;
		const ty1 = s[i + 1].type;
		const jinc = compute_jinc(ty0, ty1);
		for (let k = 0; k < jinc; k++) {
			const dk = v[j++];
			s[i].ks[k] += dk;
			norm += dk * dk;
		}
	}
	return norm;
}

function solve_spiro<C>(s: SpiroSeg<C>[], nseg: number) {
	let m: BandMat[] = [];
	let v: number[] = [];
	let perm: number[] = [];
	let nmat = count_vec(s, nseg);
	let n_alloc = nmat;
	let norm: number;

	if (nmat === 0) return 0;
	if (s[0].type !== "open" && s[0].type !== "corner") n_alloc *= 3;
	if (n_alloc < 5) n_alloc = 5;
	for (let j = 0; j < n_alloc; j++) {
		m[j] = new BandMat();
		v[j] = 0;
		perm[j] = 0;
	}
	for (let i = 0; i < 10; i++) {
		norm = spiro_iter(s, m, perm, v, nseg);
		if (norm < 1e-12) break;
	}
	return 0;
}

function findIntersection(p1: Point, c1: Point, c2: Point, p2: Point): null | Point {
	let d1: Point = { x: c1.x - p1.x, y: c1.y - p1.y };
	let d2: Point = { x: c2.x - p2.x, y: c2.y - p2.y };

	const det = d2.x * d1.y - d2.y * d1.x;
	if (Math.abs(det) < 1e-6) return null;
	const u = ((p2.y - p1.y) * d2.x - (p2.x - p1.x) * d2.y) / det;
	const v = ((p2.y - p1.y) * d1.x - (p2.x - p1.x) * d1.y) / det;
	if (u <= 0 || v <= 0) return null;
	return {
		x: p1.x + d1.x * u,
		y: p1.y + d1.y * u
	};
}

function toQuad2(
	x0: number,
	y0: number,
	x1: number,
	y1: number,
	x2: number,
	y2: number,
	x3: number,
	y3: number,
	bc: IBezierContext,
	subdivided: boolean
) {
	// construct a on-off-off-on sequence
	const sx = x0 + (3 / 4) * (x1 - x0);
	const sy = y0 + (3 / 4) * (y1 - y0);
	const fx = x3 + (3 / 4) * (x2 - x3);
	const fy = y3 + (3 / 4) * (y2 - y3);
	bc.curveTo(sx, sy, (sx + fx) / 2, (sy + fy) / 2, true);
	bc.curveTo(fx, fy, x3, y3, subdivided);
}

function toQuad4(
	x0: number,
	y0: number,
	x1: number,
	y1: number,
	x2: number,
	y2: number,
	x3: number,
	y3: number,
	bc: IBezierContext,
	subdivided: boolean
) {
	// construct a on-off-off-off-off-on sequence
	const sx = x0 + (3 / 8) * (x1 - x0);
	const sy = y0 + (3 / 8) * (y1 - y0);
	const p = (1 / 32) * (7 * x0 + 15 * x1 + 9 * x2 + x3);
	const q = (1 / 32) * (7 * y0 + 15 * y1 + 9 * y2 + y3);
	const r = (1 / 32) * (x0 + 9 * x1 + 15 * x2 + 7 * x3);
	const s = (1 / 32) * (y0 + 9 * y1 + 15 * y2 + 7 * y3);
	const fx = x3 + (3 / 8) * (x2 - x3);
	const fy = y3 + (3 / 8) * (y2 - y3);
	bc.curveTo(sx, sy, (sx + p) / 2, (sy + q) / 2, true);
	bc.curveTo(p, q, (r + p) / 2, (s + q) / 2, true);
	bc.curveTo(r, s, (r + fx) / 2, (s + fy) / 2, true);
	bc.curveTo(fx, fy, x3, y3, subdivided);
}

function spiroSegToBPath<C extends IBezierContext>(
	arc: SpiroArc,
	bc: C,
	depth: number,
	subdivided: boolean,
	isQuad: boolean,
	af: undefined | KnotCallback<C>,
	delta: number
) {
	if (arc.bend <= 1e-8) {
		bc.lineTo(arc.x1, arc.y1, subdivided);
	} else if (depth > MAX_DEPTH || arc.bend < delta) {
		const [a, b, c, d] = arc.toCubicBezier();
		if (isQuad) {
			if (arc.bend > 0.5 * delta) {
				toQuad2(a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y, bc, subdivided);
			} else {
				let pt = findIntersection(a, b, c, d);
				if (pt) {
					bc.curveTo(pt.x, pt.y, d.x, d.y, subdivided);
				} else {
					toQuad2(a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y, bc, subdivided);
				}
			}
		} else {
			bc.cubicTo(b.x, b.y, c.x, c.y, d.x, d.y, subdivided);
		}
	} else {
		const [pre, post] = arc.subdivide(1 / 2);
		spiroSegToBPath(pre, bc, depth + 1, true, isQuad, undefined, delta);
		spiroSegToBPath(post, bc, depth + 1, subdivided, isQuad, undefined, delta);
	}
	if (af) af.call(bc, arc.x0, arc.y0, arc.x1, arc.y1);
}

function run_spiro<C>(src: Knot<C>[], n: number) {
	const nseg = src[0].type === "open" ? n - 1 : n;
	const s = setup_path(src, n);
	if (nseg > 1) solve_spiro(s, nseg);
	return s;
}

function computeSegments<C>(spiros: Knot<C>[], isClosed: boolean) {
	let s: SpiroSeg<C>[];
	const n = spiros.length;
	if (!isClosed) {
		let oldStartType = spiros[0].type;
		let oldEndType = spiros[n - 1].type;
		spiros[0].type = "open";
		spiros[n - 1].type = "open_end";
		s = run_spiro(spiros, n);
		spiros[n - 1].type = oldEndType;
		spiros[0].type = oldStartType;
	} else {
		s = run_spiro(spiros, n);
	}
	return { s, n };
}
function collectSpiroArcs<C>(s: SpiroSeg<C>[], n: number, bc: C) {
	let arcs: SpiroArc[] = [];

	const nSegments = s[n - 1].type === "open_end" ? n - 1 : n;
	for (let i = 0; i < nSegments; i++) {
		const x0 = s[i].x;
		const y0 = s[i].y;
		const x1 = s[i + 1].x;
		const y1 = s[i + 1].y;
		if (i === 0) if (s[0].af) s[0].af.call(bc, x0, y0);
		arcs.push(new SpiroArc(s[i].ks, x0, y0, x1, y1));
	}
	return arcs;
}

function spiroToBPath<C extends IBezierContext>(
	s: SpiroSeg<C>[],
	n: number,
	bc: C,
	isQuad: boolean = false,
	delta: number = 1
) {
	const nSegments = s[n - 1].type === "open_end" ? n - 1 : n;
	for (let i = 0; i < nSegments; i++) {
		const x0 = s[i].x;
		const y0 = s[i].y;
		const x1 = s[i + 1].x;
		const y1 = s[i + 1].y;
		if (i === 0) {
			bc.moveTo(x0, y0);
			if (s[0].af) s[0].af.call(bc, x0, y0);
		}
		spiroSegToBPath(
			new SpiroArc(s[i].ks, x0, y0, x1, y1),
			bc,
			0,
			false,
			isQuad,
			s[i + 1].af,
			delta
		);
	}
}

export function spiroToArcs<C>(context: C, spiros: Knot<C>[], isClosed: boolean) {
	if (spiros.length < 1) return [];
	const { s, n } = computeSegments(spiros, isClosed);
	return collectSpiroArcs(s, n, context);
}

export function spiroToBezierOnContext<C extends IBezierContext>(
	spiros: Knot<C>[],
	isClosed: boolean,
	bc: C,
	isQuad: boolean = false,
	delta: number = 1
) {
	if (spiros.length < 1) return;
	const { s, n } = computeSegments(spiros, isClosed);
	spiroToBPath(s, n, bc, isQuad, delta);
}

export function spiroToBezier(spiros: Knot<IBezierContext>[], isClosed: boolean) {
	let c = new DefaultBezierContext();
	spiroToBezierOnContext(spiros, isClosed, c);
	return c.strands;
}