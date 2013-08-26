(ns cfd-clojure
  (:require [clojure.math.numeric-tower :as math]
            [clojure.core.match :as m]
            [clatrix.core :as clx]
            [clojure.core.matrix :as M])
  (:use (incanter core charts))
  (:import [org.jblas DoubleMatrix]))


(set! *warn-on-reflection* true)


;; provide discretisation linear in t and up to second order in x


(declare burgers-bc-err)
(defn discretize [f m u]
  (let [g0 (fn [u] (partition (:x-steps m) 1 u))
        h0 (fn [v] (map #(f m %) v))
        l (last u)
        g (fn [u] (conj (map #(f m %) (partition (:x-steps m) 1 u)) (first u)))
        h (fn [u] ((comp g #(concat % [l])) u))]
    (m/match [(:x-steps m) (:bc m)]
             [2 nil] (iterate g u)
             [3 nil] (iterate h u)
             [3 burgers-bc-err] (iterate (comp h0 (:bc m) g0) u)
             [3 _] (iterate (comp h0 g0 (:bc m)) u)
             ;; :else throw IllegalArgumentException
             )))


(defn discretize-2D [f m u]
  (iterate (partial f m) u))


;; Step 1: Linear Convection

(defn linear-convection [m [u_ni-1 u_ni]]
  (- u_ni (/ (* (:c m) (:dt m) (- u_ni u_ni-1)) (:dx m))))

;; Step 2: Non-Linear Convection

(defn non-linear-convection [m [u_ni-1 u_ni]]
  (- u_ni (/ (* u_ni (:dt m) (- u_ni u_ni-1)) (:dx m))))

;; Step 3: Diffusion

(defn diffusion-1D [m [u_ni-1 u_ni u_ni+1]]
  (+ u_ni (/ (* (:nu m) (:dt m) (+ u_ni+1 u_ni-1 (* -2. u_ni))) (Math/pow (:dx m) 2.))))


;; Step 4: Burgers' Equation

(defn burgers-u0 [t x nu]
  (let [f1 (+ (* -4 t) x -6.28318530717959)
        f2 (* 4. nu (+ t 1.))
        f3a (+ (* -4. t) x)
        f3 (/ (* -1. f3a f3a) f2)
        f4 (/ (* -1. f1 f1) f2)
        f5 (* -1. (+ (* -8. t) (* 2. x)) (exp f3))]
    (+ (/ (/ (* -2. nu (- f5 (* 2. f1 (exp f4))))
             f2)
          (+ (exp f4) (exp f3)))
       4.)))


;; In the given model, u0 and u101 refer to the same location modulo 2*pi.
;; Hence, the same location is used twice in updating the BC.
;; This is most likely a bug.
;;
(defn burgers-bc-err [g0u]
  (let [[u0 u1 u2] (vec (first g0u))
        [u99 u100 u101] (vec (last g0u))]
    (concat (conj g0u [u101 u0 u1]) [[u100 u101 u0]])))


;; (~ denotes 2*pi periodicy)
;; if un[0] ~ un[nx]
;; then un[1] ~ un[nx + 1]
;; and un[-1] ~ un[nx - 1]
;;
(defn burgers-bc [u]
  (concat (conj u (last (butlast u))) [(second u)]))


(defn burgers-eqn [m [u_ni-1 u_ni u_ni+1]]
  (let [dx (:dx m)
        nu (:nu m)
        dt (:dt m)]
    (+ u_ni
       (* u_ni (/ dt dx) (- u_ni-1 u_ni))
       (* nu (/ dt (* dx dx)) (+ u_ni+1 (* -2. u_ni) u_ni-1)))))


;; Step 5: 2d Linear Convection

(defn linear-convection-2D [m un]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        A (sel un :except-rows 0 :except-cols 0)
        B (sel un :except-rows upper_y :except-cols 0)
        C (sel un :except-rows 0 :except-cols upper_x)
        k (* -1. (:c m) (:dt m))
        kx (/ k (:dx m))
        ky (/ k (:dy m))
        u_core (sel (plus (mult (+ 1. kx ky) A)
                          (mult (* -1. kx) B)
                          (mult (* -1. ky) C))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        v (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set v y x (sel u_core (dec y) (dec x))))
    v))


;; Step 6: 2D Convection

;;u[1:,1:]=un[1:,1:]-(c*dt/dx*(un[1:,1:]-un[0:-1,1:]))-(c*dt/dy*(un[1:,1:]-un[1:,0:-1]))
;;u[1:,1:]=un[1:,1:]-(un[1:,1:]*dt/dx*(un[1:,1:]-un[0:-1,1:]))-vn[1:,1:]*dt/dy*(un[1:,1:]-un[1:,0:-1])
;;v[1:,1:]=vn[1:,1:]-(un[1:,1:]*dt/dx*(vn[1:,1:]-vn[0:-1,1:]))-vn[1:,1:]*dt/dy*(vn[1:,1:]-vn[1:,0:-1])

;;u[1:,1:]=Au-(c*dt/dx*(Au-Bu))-(c*dt/dy*(Au-Cu))
;;u[1:,1:]=Au-(Au*dt/dx*(Au-Bu))-Av*dt/dy*(Au-Cu)
;;v[1:,1:]=Av-(Au*dt/dx*(Av-Bv))-Av*dt/dy*(Av-Cv)

(defn convection-2D [m [un vn]]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        Au (sel un :except-rows 0 :except-cols 0)
        Bu (sel un :except-rows upper_y :except-cols 0)
        Cu (sel un :except-rows 0 :except-cols upper_x)
        Av (sel vn :except-rows 0 :except-cols 0)
        Bv (sel vn :except-rows upper_y :except-cols 0)
        Cv (sel vn :except-rows 0 :except-cols upper_x)
        kx (/ (* -1. (:dt m)) (:dx m))
        ky (/ (* -1. (:dt m)) (:dy m))
        u_core (sel (plus Au
                          (mult kx Au (minus Au Bu))
                          (mult ky Av (minus Au Cu)))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        v_core (sel (plus Av
                          (mult kx Au (minus Av Bv))
                          (mult ky Av (minus Av Cv)))
                    :except-rows (dec upper_y)
                    :except-cols (dec upper_x))
        uu (matrix 1. (:ny m) (:nx m))
        vv (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set uu y x (sel u_core (dec y) (dec x)))
      (clx/set vv y x (sel v_core (dec y) (dec x))))
    [uu vv]))


;; Step 7: 2D Diffusion

;;u[1:-1,1:-1]=un[1:-1,1:-1]
;;            +nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])
;;            +nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])

;;u[1:-1,1:-1]=D+nu*dt/dx**2*(E-2*D+F)+nu*dt/dy**2*(G-2*D+H)

;;u[1:-1,1:-1]=D+kx*(E-2*D+F)+ky*(G-2*D+H)

;;u[1:-1,1:-1] ; D :except-rows  first     last    :except-cols  first     last
;;u[2:,1:-1]   ; E :except-rows  first     second  :except-cols  first     last
;;u[0:-2,1:-1] ; F :except-rows  prev-last last    :except-cols  first     last
;;u[1:-1,0:-2] ; H :except-rows  first     last    :except-cols  prev-last last

(defn diffusion-2D [m un]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        D (sel (sel un :except-rows upper_y :except-cols upper_x )
               :except-rows 0   :except-cols 0)
        E (sel (sel un :except-rows 0 :except-cols upper_x)
               :except-rows 0 :except-cols 0)
        F (sel (sel un :except-rows upper_y :except-cols upper_x)
               :except-rows (dec upper_y) :except-cols 0)
        G (sel (sel un :except-rows upper_y :except-cols 0)
               :except-rows 0 :except-cols 0)
        H (sel (sel un :except-rows upper_y :except-cols upper_x)
               :except-rows 0 :except-cols (dec upper_x))
        k (* (:nu m) (:dt m))
        kx (/ k (math/expt (:dx m) 2.))
        ky (/ k (math/expt (:dy m) 2.))
        u_core (plus (mult (+ 1. (* -2. kx) (* -2. ky)) D)
                     (mult kx E)
                     (mult kx F)
                     (mult ky G)
                     (mult ky H))
        v (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set v y x (sel u_core (dec y) (dec x))))
    v))


;; Step 8: 2D Burgers' Equation

;; u[1:-1,1:-1] = un[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[0:-2,1:-1])-dt/dy*vn[1:-1,1:-1]* \
;;                (un[1:-1,1:-1]-un[1:-1,0:-2])+nu*dt/dx**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])+ \
;;                nu*dt/dy**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])
;;
;; v[1:-1,1:-1] = vn[1:-1,1:-1] - dt/dx*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-dt/dy*vn[1:-1,1:-1]* \
;;                (vn[1:-1,1:-1]-vn[1:-1,0:-2])+nu*dt/dx**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])+ \
;;                nu*dt/dy**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])

;; z[1:-1,1:-1] ; Dz :except-rows  first     last    :except-cols  first     last
;; z[2:,1:-1]   ; Ez :except-rows  first     second  :except-cols  first     last
;; z[0:-2,1:-1] ; Fz :except-rows  prev-last last    :except-cols  first     last
;; z[1:-1,2:]   ; Gz :except-rows  first     last    :except-cols  first     second
;; z[1:-1,0:-2] ; Hz :except-rows  first     last    :except-cols  prev-last last


;; u[1:-1,1:-1] = Du - dt/dx*Du*(Du-Fu)-dt/dy*Dv* \
;;                (Du-Hu)+nu*dt/dx**2*(Eu-2*Du+Fu)+ \
;;                nu*dt/dy**2*(Gu-2*Du+Hu)
;;
;; v[1:-1,1:-1] = Dv - dt/dx*Du*(Dv-Fv)-dt/dy*Dv* \
;;                (Dv-Hv)+nu*dt/dx**2*(Ev-2*Dv+Fv)+ \
;;                nu*dt/dy**2*(Gv-2*Dv+Hv)

(defn burgers-eqn-2D [m [un vn]]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        Du (sel (sel un :except-rows upper_y :except-cols upper_x )
               :except-rows 0   :except-cols 0)
        Eu (sel (sel un :except-rows 0 :except-cols upper_x)
               :except-rows 0 :except-cols 0)
        Fu (sel (sel un :except-rows upper_y :except-cols upper_x)
               :except-rows (dec upper_y) :except-cols 0)
        Gu (sel (sel un :except-rows upper_y :except-cols 0)
               :except-rows 0 :except-cols 0)
        Hu (sel (sel un :except-rows upper_y :except-cols upper_x)
                :except-rows 0 :except-cols (dec upper_x))
        Dv (sel (sel vn :except-rows upper_y :except-cols upper_x )
               :except-rows 0   :except-cols 0)
        Ev (sel (sel vn :except-rows 0 :except-cols upper_x)
               :except-rows 0 :except-cols 0)
        Fv (sel (sel vn :except-rows upper_y :except-cols upper_x)
               :except-rows (dec upper_y) :except-cols 0)
        Gv (sel (sel vn :except-rows upper_y :except-cols 0)
              :except-rows 0 :except-cols 0)
        Hv (sel (sel vn :except-rows upper_y :except-cols upper_x)
                :except-rows 0 :except-cols (dec upper_x))
        kx (/ (* -1. (:dt m)) (:dx m))
        kxx (/ (* (:nu m) kx) (:dx m))
        ky (/ (* -1. (:dt m)) (:dy m))
        kyy (/ (* (:nu m) ky) (:dy m))
        u_core (plus (mult (+ 1. (* 2. kxx) (* 2. kyy)) Du)
                     (mult kx Du (minus Du Fu))
                     (mult ky Dv (minus Du Hu))
                     (mult -1. kxx Eu)
                     (mult -1. kxx Fu)
                     (mult -1. kyy Gu)
                     (mult -1. kyy Hu))
        v_core (plus (mult (+ 1. (* 2. kxx) (* 2. kyy)) Dv)
                     (mult kx Du (minus Dv Fv))
                     (mult ky Dv (minus Dv Hv))
                     (mult -1. kxx Ev)
                     (mult -1. kxx Fv)
                     (mult -1. kyy Gv)
                     (mult -1. kyy Hv))
        uu (matrix 1. (:ny m) (:nx m))
        vv (matrix 1. (:ny m) (:nx m))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set uu y x (sel u_core (dec y) (dec x)))
      (clx/set vv y x (sel v_core (dec y) (dec x))))
    [uu vv]))


;;; Step 9: 2D Laplace Equation

;;  l1norm = (np.sum(np.abs(p[:])-np.abs(pn[:])))/np.sum(np.abs(pn[:]))
;;
(defn L1-norm [x y]
  (/ (reduce + (map sum (minus (abs x) (abs y)))) (reduce + (map sum (abs y)))))


;; p[2:,1:-1] ; E :except-rows first  :except-rows second :except-cols first :except-cols last
;; p[0:-2,1:-1] ; F :except-rows prev-last  :except-rows last :except-cols first :except-cols last
;; p[1:-1,2:] ; G :except-rows first  :except-rows last :except-cols first :except-cols second
;; p[1:-1,0:-2] ; H :except-rows first  :except-rows last :except-cols prev-last :except-cols last
;;
;; dxx = dx**2/(2*(dx**2+dy**2))
;; dyy = dy**2/(2*(dx**2+dy**2))
;;
;; p[1:-1,1:-1] = dyy(E+F)+dxx*(G+H)
;; p[0,0] = dyy*(pn[1,0]+pn[-1,0])+dxx*(pn[0,1]+pn[0,-1])
;; p[-1,-1] = dyy*(pn[0,-1]+pn[-2,-1])+dxx*(pn[-1,0]+pn[-1,-2])
;;
;; p[1:-1,1:-1] = (dy**2*(pn[2:,1:-1]+pn[0:-2,1:-1])+dx**2*(pn[1:-1,2:]+pn[1:-1,0:-2]))/(2*(dx**2+dy**2))
;; p[0,0] = (dy**2*(pn[1,0]+pn[-1,0])+dx**2*(pn[0,1]+pn[0,-1]))/(2*(dx**2+dy**2))
;; p[-1,-1] = (dy**2*(pn[0,-1]+pn[-2,-1])+dx**2*(pn[-1,0]+pn[-1,-2]))/(2*(dx**2+dy**2))
;;
;; p[:,0] = 0		##p = 0 @ x = 0       ;;  redundant:   p[:,0] = const.
;; p[:,-1] = y		##p = y @ x = 2       ;;  redundant:   y = pn[:,-1] = const.
;; p[0,:] = p[1,:]	##dp/dy = 0 @ y = 0
;; p[-1,:] = p[-2,:]	##dp/dy = 0 @ y = 1
;;
(defn- f-laplace [dyy dxx upper_y upper_x pn]
  (let [E (sel (sel pn :except-rows 0 :except-cols upper_x)
               :except-rows 0 :except-cols 0)
        F (sel (sel pn :except-rows upper_y :except-cols upper_x)
               :except-rows (dec upper_y) :except-cols 0)
        G (sel (sel pn :except-rows upper_y :except-cols 0)
               :except-rows 0 :except-cols 0)
        H (sel (sel pn :except-rows upper_y :except-cols upper_x)
               :except-rows 0 :except-cols (dec upper_x))
        p_core (plus (mult dyy (plus E F))
                     (mult dxx (plus G H)))
        p (matrix 0. (inc upper_y) (inc upper_x))]
    (doseq [y (range 1 upper_y)
            x (range 1 upper_x)]
      (clx/set p y x (sel p_core (dec y) (dec x))))
    (doseq [y (range (inc upper_y))]
      (clx/set p y upper_x (sel pn y upper_x)))
    (doseq [x (range (inc upper_x))]
      (clx/set p 0 x (sel p 1 x))
      (clx/set p upper_y x (sel p (dec upper_y) x)))
    p))

(defn laplace-eqn-2D [m pn]
  (let [upper_x (dec (:nx m)) ; cols
        upper_y (dec (:ny m)) ; rows
        dx2 (math/expt (:dx m) 2.)
        dy2 (math/expt (:dy m) 2.)
        edxdy (* 2. (+ dx2 dy2))
        dxx (/ dx2 edxdy)
        dyy (/ dy2 edxdy)]
    (loop [pn pn
           p (f-laplace dyy dxx upper_y upper_x pn)]
      (if (< (L1-norm p pn) (:eps m))
        p
        (recur p (f-laplace dyy dxx upper_y upper_x p))))))


;;; Step 10: 2D Poisson Equation

;; pd[2:nx,1:ny-1] ; E :except-rows first  :except-rows second :except-cols first :except-cols last
;; pd[0:nx-2,1:ny-1] ; F :except-rows prev-last  :except-rows last :except-cols first :except-cols last
;; pd[1:nx-1,2:ny] ; G :except-rows first  :except-rows last :except-cols first :except-cols second
;; pd[1:nx-1,0:ny-2] ; H :except-rows first  :except-rows last :except-cols prev-last :except-cols last
;;
;; dxx = dx**2/(2*(dx**2+dy**2))
;; dyy = dy**2/(2*(dx**2+dy**2))
;;
;; p[1:nx-1,1:ny-1] = ( dy**2/(2*(dx**2+dy**2))*(pd[2:nx,1:ny-1]+pd[0:nx-2,1:ny-1]) +
;;                      dx**2/(2*(dx**2+dy**2))*(pd[1:nx-1,2:ny]+pd[1:nx-1,0:ny-2]) -
;;                     b[1:nx-1,1:ny-1]*dx**2*dy**2/(2*(dx**2+dy**2)) )
;;
;; p[0,:] = p[nx-1,:] = p[:,0] = p[:,ny-1] = 0.0
;;
(defn poisson-eqn-2D [m pn]
  (let [upper_x (dec (:nx m)) ; rows
        upper_y (dec (:ny m)) ; cols
        E (sel (sel pn :except-rows 0 :except-cols upper_y)
               :except-rows 0 :except-cols 0)
        F (sel (sel pn :except-rows upper_x :except-cols upper_y)
               :except-rows (dec upper_x) :except-cols 0)
        G (sel (sel pn :except-rows upper_x :except-cols 0)
               :except-rows 0 :except-cols 0)
        H (sel (sel pn :except-rows upper_x :except-cols upper_y)
               :except-rows 0 :except-cols (dec upper_y))
        dx2 (math/expt (:dx m) 2.)
        dy2 (math/expt (:dy m) 2.)
        edxdy (* 2. (+ dx2 dy2))
        dxx (/ dx2 edxdy)
        dyy (/ dy2 edxdy)
        p_core (plus (mult dyy (plus E F))
                     (mult dxx (plus G H))
                     (mult (* -1. dx2 dyy) (:b m)))
        p (matrix 0. (:ny m) (:nx m))]
    (doseq [x (range 1 upper_x)
            y (range 1 upper_y)]
      (clx/set p x y (sel p_core (dec x) (dec y))))
    p))


;;; Step 11: Cavity Flow with Navier-Stokes

;; (defmacro E [z upper_x upper_y]
;;   `(sel (sel ~z :except-rows ~upper_x :except-cols ~upper_y)
;;         :except-rows 0 :except-cols 0))
;;
;; (defmacro F [z upper_x upper_y]
;;   `(sel (sel ~z :except-rows 0 :except-cols ~upper_y)
;;         :except-rows 0 :except-cols 0))
;;
;; (defmacro G [z upper_x upper_y]
;;   `(sel (sel ~z :except-rows ~upper_x :except-cols ~upper_y)
;;         :except-rows (dec ~upper_x) :except-cols 0))
;;
;; (defmacro H [z upper_x upper_y]
;;   `(sel (sel ~z :except-rows ~upper_x :except-cols 0)
;;         :except-rows 0 :except-cols 0))
;;
;; (defmacro J [z upper_x upper_y]
;;   `(sel (sel ~z :except-rows ~upper_x :except-cols ~upper_y)
;;         :except-rows 0 :except-cols (dec ~upper_y)))

(defmacro coerce [nx ny M]
   `(DoubleMatrix. ~nx ~ny ~M))


;; ;; http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
;; ;;
;; (defmacro deep-aget
;;   ([hint array idx]
;;     `(aget ~(vary-meta array assoc :tag hint) ~idx))
;;   ([hint array idx & idxs]
;;     `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~idx)]
;;        (deep-aget ~hint a# ~@idxs))))
;; ;;
;; (defmacro deep-aset [hint array & idxsv]
;;   (let [hints '{doubles double ints int} ; writing a comprehensive map is left as an exercise to the reader
;;         [v idx & sxdi] (reverse idxsv)
;;         idxs (reverse sxdi)
;;         v (if-let [h (hints hint)] (list h v) v)
;;         nested-array (if (seq idxs)
;;                        `(deep-aget ~'objects ~array ~@idxs)
;;                         array)
;;         a-sym (with-meta (gensym "a") {:tag hint})]
;;       `(let [~a-sym ~nested-array]
;;          (aset ~a-sym ~idx ~v))))
;;
;; (defmacro aget! [M i j]
;;   `(deep-aget ~'doubles ~M ~i ~j))
;;
;; (defmacro aset! [M i j v]
;;   `(deep-aset ~'doubles ~M ~i ~j ~v))

;; http://www.bestinclass.dk/index.clj/2010/03/functional-fluid-dynamics-in-clojure.html
;;
(defmacro aget!
  ([array y]      `(aget ~(vary-meta array assoc :tag 'doubles) ~y))
  ([array x & ys] `(let [a# (aget ~(vary-meta array assoc :tag 'objects) ~@ys)]
                     (aget! a# ~x))))
;;
(defmacro aset! [array x y v]
  (let [nested-array `(aget ~(vary-meta array assoc :tag 'objects) ~y)
        a-sym         (with-meta (gensym "a") {:tag 'doubles})]
    `(let [~a-sym ~nested-array]
       (aset ~a-sym ~x (double ~v)))))



(defn submatrix [M [[a nx] [b ny]]]
  (let [m (clx/matrix (DoubleMatrix. nx ny (make-array Double/TYPE (* nx ny))))]
    (dotimes [i nx]
      (dotimes [j ny]
        (clx/set m i j (clx/get M (+ a i) (+ b j)))))
    m))


(defn submatrix2 [M [[a nx] [b ny]]]
  (let [m (make-array Double/TYPE nx ny)]
    (dotimes [i nx]
      (dotimes [j ny]
        (aset! m i j (clx/get M (+ a i) (+ b j)))))
    (clx/matrix (DoubleMatrix. ^"[[D" m))))

;; (defn submatrix [M [[a len_x] [b len_y]]]
;;   (clx/from-indices len_x len_y (fn [i j] (clx/get M (+ a i) (+ b j)))))


(defmacro E [z len_x len_y]
  `(submatrix ~z [[1 ~len_x] [1 ~len_y]]))

(defmacro F [z len_x len_y]
  `(submatrix ~z [[2 ~len_x] [1 ~len_y]]))

(defmacro G [z len_x len_y]
  `(submatrix ~z [[0 ~len_x] [1 ~len_y]]))

(defmacro H [z len_x len_y]
  `(submatrix ~z [[1 ~len_x] [2 ~len_y]]))

(defmacro J [z len_x len_y]
  `(submatrix ~z [[1 ~len_x] [0 ~len_y]]))


;; (defn sc-* [x M len_x len_y]
;;   (let [xM (make-array Double/TYPE len_x len_y)]
;;     (doseq [i (range len_x)
;;             j (range len_y)]
;;       (aset! xM i j (* x (aget! M i j))))
;;     (coerce len_x len_y xM)))

(defn sc-* [x M nx ny]
  (clx/from-indices nx ny (fn [i j] (* x (clx/get M i j)))))


(defn buildup-b [m [un vn]] ; treat rho as external factor
  (let [nx-2 (- (:nx m) 2) ; rows
        ny-2 (- (:ny m) 2) ; cols
        Fu (F un nx-2 ny-2)
        Gu (G un nx-2 ny-2)
        Hu (H un nx-2 ny-2)
        Ju (J un nx-2 ny-2)
        Fv (F vn nx-2 ny-2)
        Gv (G vn nx-2 ny-2)
        Hv (H vn nx-2 ny-2)
        Jv (J vn nx-2 ny-2)
        Fu-Gu (clx/- Fu Gu)
        Hu-Ju (clx/- Hu Ju)
        Fv-Gv (clx/- Fv Gv)
        Hv-Jv (clx/- Hv Jv)]
    (clx/+ (sc-* (/  1. (* 2. (:dt m) (:dx m))) Fu-Gu nx-2 ny-2)
           (sc-* (/  1. (* 2. (:dt m) (:dy m))) Hv-Jv nx-2 ny-2)
           (sc-* (/ -1. (* 2. (:dx m) (:dx m))) (clx/* Fu-Gu Fu-Gu) nx-2 ny-2)
           (sc-* (/ -1. (* 2. (:dx m) (:dy m))) (clx/* Hu-Ju Fv-Gv) nx-2 ny-2)
           (sc-* (/ -1. (* 2. (:dy m) (:dy m))) (clx/* Hv-Jv Hv-Jv) nx-2 ny-2))))


(defn pressure-poisson [m bn pn]
  (let [nx (:nx m)   ; rows
        ny (:ny m)   ; cols
        nx-1 (dec nx)
        ny-1 (dec ny)
        nx-2 (dec nx-1)
        ny-2 (dec ny-1)
        F (F pn nx-2 ny-2)
        G (G pn nx-2 ny-2)
        H (H pn nx-2 ny-2)
        J (J pn nx-2 ny-2)
        dx2 (* (:dx m) (:dx m))
        dy2 (* (:dy m) (:dy m))
        edxdy (* 2. (+ dx2 dy2))
        dxx (/ dx2 edxdy)
        dyy (/ dy2 edxdy)
        p_core (clx/+ (sc-* dyy (clx/+ F G) nx-2 ny-2)
                      (sc-* dxx (clx/+ H J) nx-2 ny-2)
                      (sc-* (* -1. dx2 dyy) bn nx-2 ny-2))
        p (make-array Double/TYPE nx ny)]
    (doseq [x (range 1 nx-1)
            y (range 1 ny-1)]
      (clx/set p x y (clx/get p_core (dec x) (dec y))))
    (doseq [y (range ny)]
      (clx/set p 0 y (clx/get p 1 y))
      (clx/set p nx-1 y (clx/get p nx-2 y)))
    (doseq [x (range nx)]
      (clx/set p x 0 (clx/get p x 1))
      (clx/set p x ny-1 0.))
    (coerce nx ny p)))



(defn cavity-flow-2D [m [un vn pn]]
  (let [nx (:nx m)
        ny (:ny m)
        nx-1 (dec nx)
        ny-1 (dec ny)
        nx-2 (dec nx-1)
        ny-2 (dec ny-1)
        b (clx/* (:rho m) (buildup-b m [un vn]))
        ;;_ (println "before p, nt = " (:nt m))
        p (last (take (inc (:nit m)) (iterate (partial pressure-poisson m b) pn)))
        ;; p ((apply comp (repeat (inc (:nit m)) (partial pressure-poisson m b))) pn)
        ;; p (loop [n (:nit m)
        ;;          p pn]
        ;;     (if (< n 0)
        ;;       p
        ;;       (recur (dec n) ((partial pressure-poisson m b) p))))
        ;;_ (println "after p")
        Eu (E un nx-2 ny-2)
        Fu (F un nx-2 ny-2)
        Gu (G un nx-2 ny-2)
        Hu (H un nx-2 ny-2)
        Ju (J un nx-2 ny-2)
        Ev (E vn nx-2 ny-2)
        Fv (F vn nx-2 ny-2)
        Gv (G vn nx-2 ny-2)
        Hv (H vn nx-2 ny-2)
        Jv (J vn nx-2 ny-2)
        Fp (F p nx-2 ny-2)
        Gp (G p nx-2 ny-2)
        Hp (H p nx-2 ny-2)
        Jp (J p nx-2 ny-2)
        dtx (* -1. (/ (:dt m) (:dx m)))
        dty (* -1. (/ (:dt m) (:dy m)))
        dtxrho (/ dtx (* 2. (:rho m)))
        dtyrho (/ dty (* 2. (:rho m)))
        dtx2nu (/ (* (:nu m) (:dt m)) (* (:dx m) (:dx m)))
        dty2nu (/ (* (:nu m) (:dt m)) (* (:dy m) (:dy m)))
        u_core (clx/+ (sc-* (+ 1. (* -2. dtx2nu) (* -2. dty2nu)) Eu nx-2 ny-2)
                     (clx/* (sc-* dtx Eu nx-2 ny-2) (clx/- Eu Gu))
                     (clx/* dty Ev (clx/- Eu Ju))
                     (sc-* dtxrho (clx/- Fp Gp) nx-2 ny-2)
                     (sc-* dtx2nu Fu nx-2 ny-2)
                     (sc-* dtx2nu Gu nx-2 ny-2)
                     (sc-* dty2nu Hu nx-2 ny-2)
                     (sc-* dty2nu Ju nx-2 ny-2))
        v_core (clx/+ (sc-* (+ 1. (* -2. dtx2nu) (* -2. dty2nu)) Ev nx-2 ny-2)
                     (clx/* (sc-* dtx Eu nx-2 ny-2) (clx/- Ev Gv))
                     (clx/* (sc-* dty Ev nx-2 ny-2) (clx/- Ev Jv))
                     (sc-* dtyrho (clx/- Hp Jp) nx-2 ny-2)
                     (sc-* dtx2nu Fv nx-2 ny-2)
                     (sc-* dtx2nu Gv nx-2 ny-2)
                     (sc-* dty2nu Hv nx-2 ny-2)
                     (sc-* dty2nu Jv nx-2 ny-2))
        uu (make-array Double/TYPE nx ny)
        vv (make-array Double/TYPE nx ny)]
    (doseq [x (range 1 nx-1)
            y (range 1 ny-1)]
      (clx/set uu x y (clx/get u_core (dec x) (dec y)))
      (clx/set vv x y (clx/get v_core (dec x) (dec y))))
    (doseq [x (range nx-1)]
      (clx/set uu x ny-1 1.))
    [(coerce nx ny uu)
     (coerce nx ny vv)
     (coerce nx ny p)]))
