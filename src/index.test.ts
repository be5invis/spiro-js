import test from "ava";

import { spiroToArcs } from "./index";

test("Derivative and eval works", t => {
	const arcs = spiroToArcs(
		null,
		[
			{ type: "g4", x: 0, y: 0 },
			{ type: "g4", x: 100, y: 0 },
			{ type: "g4", x: 100, y: 100 },
			{ type: "g4", x: 0, y: 100 }
		],
		true
	);

	t.true(arcs.length > 0);

	const arc = arcs[0];
	const stops = 16;
	for (let stop = 0; stop <= stops; stop++) {
		const s = stop / stops;
		const [pre] = arc.subdivide(s);
		const z = arc.eval(s);
		const dz = arc.derivative(s);

		// End positions are the same
		t.true(Math.hypot(z.x - pre.x1, z.y - pre.y1) < 1e-8);

		// Derivatives are proportional to the arc length
		t.true(
			Math.hypot(
				(dz.x * pre.arcLength) / arc.arcLength - pre.deriveX1,
				(dz.y * pre.arcLength) / arc.arcLength - pre.deriveY1
			) < 1e-8
		);
	}
});
