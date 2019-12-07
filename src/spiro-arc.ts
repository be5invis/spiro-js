import { Point } from "./context";

const N = 4;
export type SpiroK = [number, number, number, number];

// Integrate polynomial spiral curve over range -.5 .. .5.
export function integrate_spiro(ks: SpiroK): [number, number] {
	const th1 = ks[0];
	const th2 = 0.5 * ks[1];
	const th3 = (1 / 6) * ks[2];
	const th4 = (1 / 24) * ks[3];
	const ds = 1 / N;
	const ds2 = ds * ds;
	const ds3 = ds2 * ds;
	const k0 = ks[0] * ds;
	const k1 = ks[1] * ds;
	const k2 = ks[2] * ds;
	const k3 = ks[3] * ds;

	let s = 0.5 * ds - 0.5;

	let x = 0;
	let y = 0;

	for (let i = 0; i < N; i++) {
		let u: number, v: number;
		let km0: number, km1: number, km2: number, km3: number;

		km0 = (((1 / 6) * k3 * s + 0.5 * k2) * s + k1) * s + k0;
		km1 = ((0.5 * k3 * s + k2) * s + k1) * ds;
		km2 = (k3 * s + k2) * ds2;
		km3 = k3 * ds3;

		// Order 12 implementation
		const t1_1 = km0;
		const t1_2 = 0.5 * km1;
		const t1_3 = (1 / 6) * km2;
		const t1_4 = (1 / 24) * km3;
		const t2_2 = t1_1 * t1_1;
		const t2_3 = 2 * (t1_1 * t1_2);
		const t2_4 = 2 * (t1_1 * t1_3) + t1_2 * t1_2;
		const t2_5 = 2 * (t1_1 * t1_4 + t1_2 * t1_3);
		const t2_6 = 2 * (t1_2 * t1_4) + t1_3 * t1_3;
		const t2_7 = 2 * (t1_3 * t1_4);
		const t2_8 = t1_4 * t1_4;
		const t3_4 = t2_2 * t1_2 + t2_3 * t1_1;
		const t3_6 = t2_2 * t1_4 + t2_3 * t1_3 + t2_4 * t1_2 + t2_5 * t1_1;
		const t3_8 = t2_4 * t1_4 + t2_5 * t1_3 + t2_6 * t1_2 + t2_7 * t1_1;
		const t3_10 = t2_6 * t1_4 + t2_7 * t1_3 + t2_8 * t1_2;
		const t4_4 = t2_2 * t2_2;
		const t4_5 = 2 * (t2_2 * t2_3);
		const t4_6 = 2 * (t2_2 * t2_4) + t2_3 * t2_3;
		const t4_7 = 2 * (t2_2 * t2_5 + t2_3 * t2_4);
		const t4_8 = 2 * (t2_2 * t2_6 + t2_3 * t2_5) + t2_4 * t2_4;
		const t4_9 = 2 * (t2_2 * t2_7 + t2_3 * t2_6 + t2_4 * t2_5);
		const t4_10 = 2 * (t2_2 * t2_8 + t2_3 * t2_7 + t2_4 * t2_6) + t2_5 * t2_5;
		const t5_6 = t4_4 * t1_2 + t4_5 * t1_1;
		const t5_8 = t4_4 * t1_4 + t4_5 * t1_3 + t4_6 * t1_2 + t4_7 * t1_1;
		const t5_10 = t4_6 * t1_4 + t4_7 * t1_3 + t4_8 * t1_2 + t4_9 * t1_1;
		const t6_6 = t4_4 * t2_2;
		const t6_7 = t4_4 * t2_3 + t4_5 * t2_2;
		const t6_8 = t4_4 * t2_4 + t4_5 * t2_3 + t4_6 * t2_2;
		const t6_9 = t4_4 * t2_5 + t4_5 * t2_4 + t4_6 * t2_3 + t4_7 * t2_2;
		const t6_10 = t4_4 * t2_6 + t4_5 * t2_5 + t4_6 * t2_4 + t4_7 * t2_3 + t4_8 * t2_2;
		const t7_8 = t6_6 * t1_2 + t6_7 * t1_1;
		const t7_10 = t6_6 * t1_4 + t6_7 * t1_3 + t6_8 * t1_2 + t6_9 * t1_1;
		const t8_8 = t6_6 * t2_2;
		const t8_9 = t6_6 * t2_3 + t6_7 * t2_2;
		const t8_10 = t6_6 * t2_4 + t6_7 * t2_3 + t6_8 * t2_2;
		const t9_10 = t8_8 * t1_2 + t8_9 * t1_1;
		const t10_10 = t8_8 * t2_2;
		u = 1;
		v = 0;
		v += (1 / 12) * t1_2 + (1 / 80) * t1_4;
		u -= (1 / 24) * t2_2 + (1 / 160) * t2_4 + (1 / 896) * t2_6 + (1 / 4608) * t2_8;
		v -= (1 / 480) * t3_4 + (1 / 2688) * t3_6 + (1 / 13824) * t3_8 + (1 / 67584) * t3_10;
		u += (1 / 1920) * t4_4 + (1 / 10752) * t4_6 + (1 / 55296) * t4_8 + (1 / 270336) * t4_10;
		v += (1 / 53760) * t5_6 + (1 / 276480) * t5_8 + (1 / 1.35168e6) * t5_10;
		u -= (1 / 322560) * t6_6 + (1 / 1.65888e6) * t6_8 + (1 / 8.11008e6) * t6_10;
		v -= (1 / 1.16122e7) * t7_8 + (1 / 5.67706e7) * t7_10;
		u += (1 / 9.28973e7) * t8_8 + (1 / 4.54164e8) * t8_10;
		v += (1 / 4.08748e9) * t9_10;
		u -= (1 / 4.08748e10) * t10_10;

		const th = (((th4 * s + th3) * s + th2) * s + th1) * s;
		const cth = Math.cos(th);
		const sth = Math.sin(th);

		x += cth * u - sth * v;
		y += cth * v + sth * u;
		s += ds;
	}

	return [x * ds, y * ds];
}

export function ksBend(ks: SpiroK) {
	return (
		Math.abs(ks[0]) +
		Math.abs(0.5 * ks[1]) +
		Math.abs(0.125 * ks[2]) +
		Math.abs((1 / 48) * ks[3])
	);
}

export function ksTheta(s0: number, s1: number, ks: SpiroK) {
	const s = (s0 + s1) / 2;
	return (
		s * ks[0] +
		(1 / 2) * s * s * ks[1] +
		(1 / 6) * s * s * s * ks[2] +
		(1 / 24) * s * s * s * s * ks[3]
	);
}

function divideKs(s0: number, s1: number, ks: SpiroK): SpiroK {
	const s = (s0 + s1) / 2;
	const t = s1 - s0;
	return [
		t * (ks[0] + s * ks[1] + (1 / 2) * s * s * ks[2] + (1 / 6) * s * s * s * ks[3]),
		t * t * (ks[1] + s * ks[2] + (1 / 2) * s * s * ks[3]),
		t * t * t * (ks[2] + s * ks[3]),
		t * t * t * t * ks[3]
	];
}

// A spiro arc is defined as
// ∫[-1/2 .. 1/2] exp(i * (k0 s + (1/2) * k1 * s^2 + (1/6) * k2 * s^3 + (1/24) * k3 * s^4)) d s
export class SpiroArc {
	private rot: number;

	public deriveX0: number;
	public deriveY0: number;
	public deriveX1: number;
	public deriveY1: number;

	public arcLength: number;
	public bend: number;

	constructor(
		private readonly ks: SpiroK,
		public readonly x0: number,
		public readonly y0: number,
		public readonly x1: number,
		public readonly y1: number
	) {
		const seg_ch = Math.hypot(x1 - x0, y1 - y0);
		const seg_th = Math.atan2(y1 - y0, x1 - x0);

		const xy = integrate_spiro(ks);
		const ch = Math.hypot(xy[0], xy[1]);
		const th = Math.atan2(xy[1], xy[0]);
		this.arcLength = seg_ch / ch;
		this.rot = seg_th - th;

		const thetaLeft = this.rot + ksTheta(-1 / 2, -1 / 2, ks);
		const thetaRight = this.rot + ksTheta(1 / 2, 1 / 2, ks);
		this.deriveX0 = this.arcLength * Math.cos(thetaLeft);
		this.deriveY0 = this.arcLength * Math.sin(thetaLeft);
		this.deriveX1 = this.arcLength * Math.cos(thetaRight);
		this.deriveY1 = this.arcLength * Math.sin(thetaRight);

		this.bend = ksBend(ks);
	}

	public toCubicBezier(): [Point, Point, Point, Point] {
		return [
			{ x: this.x0, y: this.y0 },
			{ x: this.x0 + this.deriveX0 / 3, y: this.y0 + this.deriveY0 / 3 },
			{ x: this.x1 - this.deriveX1 / 3, y: this.y1 - this.deriveY1 / 3 },
			{ x: this.x1, y: this.y1 }
		];
	}

	private evalImpl(t: number, kSub: SpiroK): Point {
		const thSub = this.rot + ksTheta(-1 / 2, t, this.ks);
		const cth = (t + 1 / 2) * this.arcLength * Math.cos(thSub);
		const sth = (t + 1 / 2) * this.arcLength * Math.sin(thSub);
		const xySub = integrate_spiro(kSub);
		const xMid = this.x0 + cth * xySub[0] - sth * xySub[1];
		const yMid = this.y0 + cth * xySub[1] + sth * xySub[0];
		return { x: xMid, y: yMid };
	}

	public eval(at: number): Point {
		const t = at - 1 / 2;
		const kSub = divideKs(-1 / 2, t, this.ks);
		return this.evalImpl(t, kSub);
	}

	public derivative(at: number) {
		const t = at - 1 / 2;
		const theta = this.rot + ksTheta(t, t, this.ks);
		return { x: this.arcLength * Math.cos(theta), y: this.arcLength * Math.sin(theta) };
	}

	public subdivide(at: number): [SpiroArc, SpiroArc] {
		const t = at - 1 / 2;
		const kSub = divideKs(-1 / 2, t, this.ks);
		const kSubRear = divideKs(t, 1 / 2, this.ks);

		const mid = this.evalImpl(t, kSub);

		return [
			new SpiroArc(kSub, this.x0, this.y0, mid.x, mid.y),
			new SpiroArc(kSubRear, mid.x, mid.y, this.x1, this.y1)
		];
	}
}
