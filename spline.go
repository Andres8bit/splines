package spline

import "math"

type Num float64
type Spline struct {
	z []Num
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

//Constructs the zi matrix of a cubic spline and stores the coed in z
func (s *Spline) New(n int, x Num, t, y []Num) {
	s.z = coef(n, t, y)
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
