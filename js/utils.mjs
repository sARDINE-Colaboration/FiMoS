export default class Vector extends Array {
    add(other) {
      return this.map((e, i) => e + other[i]);
    }
    mul(other) {
      return this.map((e, i) => e * other[i]);
    }
    div(other) {
      return this.map((e, i) => e / other[i]);
    }
    dot(other) {
        return this.reduce((total, current, i) => total + current * other[i], 0);
    }
    add_scalar(scalar) {
      return this.map(e => e + scalar);
    }
    mul_scalar(scalar) {
      return this.map(e => e * scalar);
    }
    div_scalar(scalar) {
      return this.map(e => e / scalar);
    }
    sum() {
        return this.reduce((total, current) => total + current, 0);
    }
    norm() {
        const sumOfSquares = this.reduce((total, current) => total + current ** 2, 0);
        return Math.sqrt(sumOfSquares);
    }
    normXY() {
        return Math.sqrt(this[0] ** 2 + this[1] ** 2);
    }
    phi() {
        return Math.atan2(this[1], this[0]);
    }
    theta() {
        return Math.atan2(this.normXY(), this[2]);
    }
    u_v() {
        return this.div_scalar(this.norm());
    }
    u_phi() {
        const phi = this.phi();
        const sin_theta = Math.sin(this.theta());
        return [-Math.sin(phi) * sin_theta, Math.cos(phi) * sin_theta, 0];
    }
    u_theta() {
        const phi = this.phi();
        const theta = this.theta();
        return [Math.cos(phi) * Math.cos(theta), Math.sin(phi) * Math.cos(theta), - Math.cos(theta)];
    }
  }