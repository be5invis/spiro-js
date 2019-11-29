export interface IBezierContext {
	moveTo(x: number, y: number): void;
	lineTo(x: number, y: number, subdivided: boolean): void;
	curveTo(x1: number, y1: number, x2: number, y2: number, subdivided: boolean): void;
	cubicTo(
		x1: number,
		y1: number,
		x2: number,
		y2: number,
		x: number,
		y: number,
		subdivided: boolean
	): void;
}
export type KnotCallback<C> = (this: C, x0: number, y0: number, x1?: number, y1?: number) => void;

export type Point = { x: number; y: number };

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

	moveTo(x: number, y: number) {
		this.lastx = x;
		this.lasty = y;
	}

	lineTo(x: number, y: number) {
		const arc: LineArc = {
			order: 1,
			start: { x: this.lastx, y: this.lasty },
			end: { x: x, y: y }
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
			end: { x: x, y: y }
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
			end: { x: x, y: y }
		};
		this.strands.push(arc);
		this.lastx = x;
		this.lasty = y;
	}
}
