export class Fish {
  constructor(id) {
    this.id = id;
    this.state = 0;
    this.positions = [];
    this.timestamp = [];
    this.timestamp_states = [];
    this.parameters = [];
  }

  initialize(position, time = 0) {
    this.positions = [position];
    this.timestamp = [time];
    this.timestamp_states = [time];
  }

  setState(new_state) {
    this.state = new_state;
  }

  step() {
    // TODO: integrate with physics and environment forces.
  }
}
