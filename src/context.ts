import { Point } from "./base";
import { SpiroArc } from "./spiro-arc";

// Arc context
export interface IArcContext {
	beginShape(): void;
	endShape(): void;
	moveTo(x: number, y: number): void;
	arcTo(arc: SpiroArc, x: number, y: number, isStraight: boolean): void;
}
export class SimplyCollectArcContext implements IArcContext {
	public arcs: SpiroArc[] = [];
	beginShape() {}
	endShape() {}
	moveTo() {}
	arcTo(arc: SpiroArc) {
		this.arcs.push(arc);
	}
}

// Bezier context
export interface IBezierContext {
	beginShape(): void;
	endShape(): void;
	moveTo(x: number, y: number): void;
	lineTo(x: number, y: number): void;
	curveTo(x1: number, y1: number, x2: number, y2: number): void;
	cubicTo(x1: number, y1: number, x2: number, y2: number, x: number, y: number): void;
}
export type KnotCallback<C> = (this: C, x0: number, y0: number, x1?: number, y1?: number) => void;

export type LineArc = {
	order: 1;
	start: Point;
	end: Point;
};
export type QuadArc = {
	order: 2;
	start: Point;
	c1: Point;
	end: Point;
};
export type CubeArc = {
	order: 3;
	start: Point;
	c1: Point;
	c2: Point;
	end: Point;
};
export type Arc = LineArc | QuadArc | CubeArc;
export class DefaultBezierContext implements IBezierContext {
	strands: Arc[] = [];
	lastx = 0;
	lasty = 0;

	beginShape() {}
	endShape() {}

	moveTo(x: number, y: number) {
		this.lastx = x;
		this.lasty = y;
	}

	lineTo(x: number, y: number) {
		const arc: LineArc = {
			order: 1,
			start: { x: this.lastx, y: this.lasty },
			end: { x: x, y: y },
		};
		this.strands.push(arc);
		this.lastx = x;
		this.lasty = y;
	}
	curveTo(x1: number, y1: number, x: number, y: number) {
		const arc: QuadArc = {
			order: 2,
			start: { x: this.lastx, y: this.lasty },
			c1: { x: x1, y: y1 },
			end: { x: x, y: y },
		};
		this.strands.push(arc);
		this.lastx = x;
		this.lasty = y;
	}
	cubicTo(x1: number, y1: number, x2: number, y2: number, x: number, y: number) {
		const arc: CubeArc = {
			order: 3,
			start: { x: this.lastx, y: this.lasty },
			c1: { x: x1, y: y1 },
			c2: { x: x2, y: y2 },
			end: { x: x, y: y },
		};
		this.strands.push(arc);
		this.lastx = x;
		this.lasty = y;
	}
}
