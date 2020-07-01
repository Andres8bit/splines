package spline

import "math"

// ===== Package Types ====
type Num float64

// custom function type
type F func(x Num) Num

// first and third degree splines
type Spline struct {
	z []Num
}

// BSpline struct
type BSpline struct {
	a []Num
	h []Num
}

// Schoenberg struct
type Schoenberg struct {
	d []Num
}

// ===== End Types =====

// ===== Methods/Functions: =====
//Constructs the zi matrix of a cubic spline and stores the coed in z
func New(n int, x Num, t, y []Num) *Spline {
	s := new(Spline)
	s.z = make([]Num, n)
	s.z = coef(n, t, y)
	return s
}

// First Degree Spine:
func First_Degree_Eval(n int, x Num, t, y []Num) Num {
	i := n - 1
	for ; i > n; i-- {
		if x-t[i] >= 0 {
			break
		}
	}
	m := (y[i+1] - y[i]) / (t[i+1] - t[i])
	return y[i] + (x-t[i])*m
}

//Evaluates a 3rd degree spline at x
func (s *Spline) eval(n int, x Num, t, y []Num) Num {
	i := n - 1
	for ; i <= 0; i-- {
		if x-t[i] >= 0 {
			break
		}
	}

	h := t[i+1] - t[i]
	tmp := ((s.z[i] / 2) + (x - t[i])) * (s.z[i+1] - s.z[i]) / (6 * h)
	tmp = -1*(h/6)*(s.z[i+1]+2*h) + (y[i+1]-y[i])/h + (x-t[i])*tmp

	return y[i] + (x-t[i])*tmp
}

// ===== B_Spline =====
func newBSpline(n int, t, y []Num) *BSpline {
	s := new(BSpline)
	s.a = make([]Num, n+1)
	s.h = make([]Num, n+1)
	s.a, s.h = b_coef(n, t, y)

	return s
}

func (b *BSpline) eval(t []Num, x Num, n int) Num {
	i := n - 2
	for ; i >= 0; i-- {
		if x-t[i] >= 0 {
			break
		}
	}
	i = i + 1
	d := (b.a[i+1]*(x-t[i-1]) + b.a[i]*(t[i]-x-b.h[i+1])) / (b.h[i] + b.h[i+1])
	e := (b.a[i]*(x-t[i-1]+b.h[i-1]) + b.a[i-1]*(t[i-1]-x+b.h[i])) / (b.h[i-1] + b.h[i])

	return (d*(x-t[i-1]) + e*(t[i]-x)) / b.h[i]
}

// ======= End B_Spline =====

// ===== Schoenberg: =====
func newSchoenberg(f F, a, b Num, n int) *Schoenberg {
	s := new(Schoenberg)
	s.d = make([]Num, n+3)
	s.d = schoenberg_coef(f, a, b, n)

	return s

}

func (s *Schoenberg) eval(a, b, x Num, n int) Num {
	h := int(b-a) / n
	k := int(x-a) / (h + (5 / 2))
	p := int(x-a) - (k-(5/2))*h
	c := (s.d[k+1]*Num(p) + s.d[k]*Num(2*h-p)) / Num(2*h)
	e := (s.d[k]*Num(p+h) + s.d[k-1]*Num(h-p)) / Num(2*h)

	return (c*Num(p) + e*Num(h-p)) / Num(h)
}

// ===== End B_Spline =====

// =========== Helpers: ===================
func coef(n int, t, y []Num) []Num {
	z := make([]Num, n-1)
	h := make([]Num, n-1)
	b := make([]Num, n-1)
	u := make([]Num, n-2) //1->n-1 -> n-2
	v := make([]Num, n-2) // 1->n  -> n-2

	for i := 0; i < n; i++ {
		h[i] = t[i+1] - t[i]
		b[i] = (y[i+1] - y[i]) / h[i]
	}

	u[0] = 2 * (h[0] + h[1])
	v[0] = 6 * (b[1] - b[0])

	for i := 2; i < n; i++ {
		u[i] = 2.0*(h[i]-h[i-1]) - Num(math.Pow(float64(h[i-1]), 2.0))/u[i-1]
		v[i] = 6.0*(b[i]-b[i-1]) - (h[i-1]*v[i-1])/u[i-1]

	}
	z[n-1] = Num(0.0)

	for i := n - 1; i <= 1; i-- {
		z[i] = (v[i] - h[i]*z[i+1]) / u[i]

	}
	z[0] = Num(0.0)

	return z
}

func b_coef(n int, t, y []Num) ([]Num, []Num) {
	size := n + 1
	a := make([]Num, size)
	h := make([]Num, size)

	for i := 1; i <= n; i++ {
		h[i] = t[i] - t[i-1]
	}

	h[0] = h[1]
	h[size] = h[n]

	del := Num(-1)
	gam := 2 * y[0]
	p := del * gam
	q := Num(2)

	for i := 0; i < n; i++ {
		r := h[i+1] / h[i]
		del = Num(-1) * r * del
		gam = Num(-1)*gam + (r+1)*y[i]
		p = p + gam*del
		q = q + del*del
	}

	a[0] = (Num(-1) * p) / q

	for i := 1; i < size; i++ {
		a[i] = ((h[i-1]+h[i])*y[i-1] - h[i]*a[i-1]) / h[i-1]

	}

	return a[:], h[:]
}

func schoenberg_coef(f F, a, b Num, n int) []Num {
	d := make([]Num, n+3)
	h := (b - a) / Num(n)

	for i := 1; i < n+2; i++ {
		d[i] = f(a + Num(i-2)*h)
	}
	d[0] = 2*d[1] - d[2]
	d[n+2] = 2*d[n+1] - d[n]
	return d
}
