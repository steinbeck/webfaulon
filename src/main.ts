import Alpine from 'alpinejs';
import { appComponent } from './ui/app';

// Register Alpine.js components
Alpine.data('app', appComponent);

// Start Alpine
Alpine.start();

// Log initialization for verification
console.log('WebFaulon loaded - Alpine.js initialized');
